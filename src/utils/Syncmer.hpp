#pragma once

#include <array>
#include <cstdint>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "syncmer_hash.hpp"

namespace taxicf::syncmer
{
namespace detail
{
inline std::vector<int> convert_positions(std::span<const size_t> positions,
                                         size_t const smer_size,
                                         size_t const kmer_size)
{
    if (smer_size == 0)
        throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
    if (kmer_size <= smer_size)
        throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};
    if (positions.empty())
        throw std::invalid_argument{"syncmer positions must not be empty"};

    std::vector<int> converted_positions;
    converted_positions.reserve(positions.size());
    for (size_t const pos : positions)
    {
        if (pos >= kmer_size - smer_size + 1)
            throw std::invalid_argument{"syncmer position out of window range"};
        converted_positions.push_back(static_cast<int>(pos));
    }

    return converted_positions;
}

template <std::ranges::forward_range rng_t>
inline void compute_hashes_append_seqan(rng_t const & sequence,
                                        size_t const smer_size,
                                        size_t const kmer_size,
                                        std::span<const size_t> positions,
                                        uint64_t const seed,
                                        bool const canonical,
                                        std::vector<uint64_t> & out)
{
    auto const converted_positions = convert_positions(positions, smer_size, kmer_size);

    if (canonical)
    {
        auto view = seqan3::detail::syncmer_hash_fn{}(sequence,
                                                      smer_size,
                                                      kmer_size,
                                                      converted_positions,
                                                      seqan3::seed{seed});
        for (auto const value : view)
            out.push_back(static_cast<uint64_t>(value));
    }
    else
    {
        auto view = seqan3::detail::syncmer_hash_no_reverse_fn{}(sequence,
                                                                 smer_size,
                                                                 kmer_size,
                                                                 converted_positions);
        for (auto const value : view)
            out.push_back(static_cast<uint64_t>(value));
    }
}

struct min_queue
{
    // kept for ABI stability (unused)
    static constexpr size_t capacity = 64;
};

struct sliding_window_min
{
    static constexpr size_t capacity = 64; // must be power-of-two
    static constexpr size_t mask = capacity - 1;

    std::array<uint64_t, capacity> values{};
    size_t head{0};      // index of the oldest element
    size_t filled{0};    // number of elements currently stored (<= window_size)
    size_t window_size{0};
    size_t min_offset{0}; // offset from head [0..window_size-1]
    uint64_t min_value{0};
    bool keep_latest_on_tie{false};

    void reset(size_t const w, bool const keep_latest) noexcept
    {
        head = 0;
        filled = 0;
        window_size = w;
        min_offset = 0;
        min_value = 0;
        keep_latest_on_tie = keep_latest;
    }

    void recompute_min() noexcept
    {
        min_offset = 0;
        min_value = values[head];
        for (size_t off = 1; off < window_size; ++off)
        {
            uint64_t const v = values[(head + off) & mask];
            if (keep_latest_on_tie ? (v <= min_value) : (v < min_value))
            {
                min_value = v;
                min_offset = off;
            }
        }
    }

    void push(uint64_t const v) noexcept
    {
        if (filled < window_size)
        {
            size_t const idx = (head + filled) & mask;
            values[idx] = v;

            if (filled == 0)
            {
                min_value = v;
                min_offset = 0;
            }
            else if (keep_latest_on_tie ? (v <= min_value) : (v < min_value))
            {
                min_value = v;
                min_offset = filled;
            }

            ++filled;
            return;
        }

        // Slide window by one: drop oldest at head, append new value at tail.
        head = (head + 1) & mask;

        bool const dropped_min = (min_offset == 0);
        if (!dropped_min)
            --min_offset;

        size_t const tail = (head + window_size - 1) & mask;
        values[tail] = v;

        if (dropped_min)
        {
            recompute_min();
        }
        else if (keep_latest_on_tie ? (v <= min_value) : (v < min_value))
        {
            min_value = v;
            min_offset = window_size - 1;
        }
    }
};

template <std::ranges::forward_range rng_t>
inline void compute_hashes_append_dna4_fast(rng_t const & sequence,
                                            size_t const smer_size,
                                            size_t const kmer_size,
                                            std::span<const size_t> positions,
                                            uint64_t const seed,
                                            std::vector<uint64_t> & out)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(std::same_as<alphabet_t, seqan3::dna4>, "dna4 fast path requires seqan3::dna4");
    static_assert(std::ranges::sized_range<rng_t>, "dna4 fast path requires sized range");

    if (smer_size == 0)
        throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
    if (kmer_size <= smer_size)
        throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};
    if (positions.empty())
        throw std::invalid_argument{"syncmer positions must not be empty"};

    size_t const window_size = kmer_size - smer_size + 1;
    if (window_size > 63)
        throw std::invalid_argument{"dna4 fast path supports window_size <= 63"};

    uint64_t pos_mask = 0;
    for (size_t const pos : positions)
    {
        if (pos >= window_size)
            throw std::invalid_argument{"syncmer position out of window range"};
        pos_mask |= (1ULL << pos);
    }

    size_t const n = std::ranges::size(sequence);
    if (n < kmer_size)
        throw std::invalid_argument{"The given sequence is too short to satisfy the given window_size.\n"
                                    "Please choose a smaller window_size."};

    if (kmer_size > 32 || smer_size > 32)
        throw std::invalid_argument{"dna4 fast path supports k<=32 and s<=32"};

    out.reserve(out.size() + (n / window_size) + 16);

    auto mask_bits = [](size_t const mer_size) -> uint64_t {
        if (mer_size == 32)
            return ~uint64_t{0};
        return (uint64_t{1} << (2 * mer_size)) - 1;
    };

    uint64_t const k_mask = mask_bits(kmer_size);
    uint64_t const s_mask = mask_bits(smer_size);
    uint64_t const k_high_shift = 2 * (kmer_size - 1);
    uint64_t const s_high_shift = 2 * (smer_size - 1);

    sliding_window_min f_min{};
    f_min.reset(window_size, false);
    sliding_window_min r_min{};
    r_min.reset(window_size, true);

    uint64_t f_kmer = 0;
    uint64_t r_kmer = 0;
    uint64_t f_smer = 0;
    uint64_t r_smer = 0;

    auto update_mer_states = [&](uint8_t const r, uint8_t const rc) noexcept {
        f_smer = ((f_smer << 2) | r) & s_mask;
        r_smer = (r_smer >> 2) | (uint64_t{rc} << s_high_shift);
        r_smer &= s_mask;

        f_kmer = ((f_kmer << 2) | r) & k_mask;
        r_kmer = (r_kmer >> 2) | (uint64_t{rc} << k_high_shift);
        r_kmer &= k_mask;
    };

    // 1) Warm up until the first s-mer.
    for (size_t pos = 0; pos + 1 < smer_size; ++pos)
    {
        uint8_t const r = static_cast<uint8_t>(seqan3::to_rank(sequence[pos]));
        update_mer_states(r, static_cast<uint8_t>(3u - r));
    }

    // 2) Build up the s-mer window until the first k-mer is available.
    for (size_t pos = smer_size - 1; pos + 1 < kmer_size; ++pos)
    {
        uint8_t const r = static_cast<uint8_t>(seqan3::to_rank(sequence[pos]));
        update_mer_states(r, static_cast<uint8_t>(3u - r));

        f_min.push(f_smer ^ seed);
        r_min.push(r_smer ^ seed);
    }

    // 3) Main loop: for each new base, update s-mer minima and emit k-mer syncmers.
    for (size_t pos = kmer_size - 1; pos < n; ++pos)
    {
        uint8_t const r = static_cast<uint8_t>(seqan3::to_rank(sequence[pos]));
        update_mer_states(r, static_cast<uint8_t>(3u - r));

        f_min.push(f_smer ^ seed);
        r_min.push(r_smer ^ seed);

        size_t const f_offset = f_min.min_offset;
        size_t const r_offset_rev = (window_size - 1) - r_min.min_offset;

        uint64_t const f_hash = f_kmer ^ seed;
        uint64_t const r_hash = r_kmer ^ seed;

        if (f_hash < r_hash)
        {
            if (pos_mask & (1ULL << f_offset))
                out.push_back(f_hash);
        }
        else
        {
            if (pos_mask & (1ULL << r_offset_rev))
                out.push_back(r_hash);
        }
    }
}

} // namespace detail

struct SyncmerOptions
{
    size_t kmer_size{};
    size_t smer_size{};
    std::vector<size_t> positions{};
    uint64_t seed{0x8F3F73B5CF1C9ADEULL};
    bool canonical{true};
};

template <std::ranges::forward_range rng_t>
inline void compute_hashes_append(rng_t const &sequence,
                                  size_t const smer_size,
                                  size_t const kmer_size,
                                  std::span<const size_t> positions,
                                  uint64_t const seed,
                                  bool const canonical,
                                  std::vector<uint64_t> &out);

template <std::ranges::forward_range rng_t>
inline std::vector<uint64_t> compute_hashes(rng_t const & sequence,
                                            size_t const smer_size,
                                            size_t const kmer_size,
                                            std::span<const size_t> positions,
                                            uint64_t const seed = 0x8F3F73B5CF1C9ADEULL,
                                            bool const canonical = true)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

    std::vector<uint64_t> values;
    compute_hashes_append(sequence, smer_size, kmer_size, positions, seed, canonical, values);
    return values;
}

inline std::vector<uint64_t> compute_hashes(seqan3::dna5_vector const & sequence,
                                            SyncmerOptions const & options)
{
    return compute_hashes(sequence,
                          options.smer_size,
                          options.kmer_size,
                          options.positions,
                          options.seed,
                          options.canonical);
}

template <std::ranges::forward_range rng_t>
inline void compute_hashes_append(rng_t const &sequence,
                                  size_t const smer_size,
                                  size_t const kmer_size,
                                  std::span<const size_t> positions,
                                  uint64_t const seed,
                                  bool const canonical,
                                  std::vector<uint64_t> &out)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

    if constexpr (std::same_as<alphabet_t, seqan3::dna4>)
    {
        if (canonical && kmer_size <= 32 && smer_size <= 32)
        {
            detail::compute_hashes_append_dna4_fast(sequence,
                                                    smer_size,
                                                    kmer_size,
                                                    positions,
                                                    seed,
                                                    out);
            return;
        }
    }

    detail::compute_hashes_append_seqan(sequence,
                                        smer_size,
                                        kmer_size,
                                        positions,
                                        seed,
                                        canonical,
                                        out);
}

} // namespace taxicf::syncmer
