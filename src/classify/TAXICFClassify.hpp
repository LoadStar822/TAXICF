/*
 * -----------------------------------------------------------------------------
 * Filename:      TAXICFClassify.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-08-10
 *
 * Last Modified: 2024-11-18
 *
 * Description:
 *  This is TAXICFClassify module for TAXICF
 *
 * Version:
 *  1.3
 * -----------------------------------------------------------------------------
 */
#pragma once

#include "classifyConfig.hpp"

#include <functional>
#include <unordered_map>
#include <vector>

namespace TAXICFClassify {
void run(ClassifyConfig config);
void postEmDecision(
    std::vector<classifyResult> &results, const DecisionConfig &decisionConfig,
    const std::unordered_map<std::string, double> &classWeights);
} // namespace TAXICFClassify
