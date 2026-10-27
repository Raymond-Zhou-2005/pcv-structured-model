# 主文与仓库迁移覆盖计划

基准：`Teacher_Revision/main.tex`（324行、27页）与 `appendix.tex`（1670行、48页）。以下行号均指基准，便于最少量修改；不是新稿的最终行号。

## 必须消除或改写的SI依赖

| 原main行 | 原依赖/说明 | 最小改法 |
|---:|---|---|
| 43 | earlier comparisons documented in supporting information | 改为仓库 `docs/extended_research_record.md` 的历史比较；早期队列本身不是外部发表自引 |
| 108 | Supplementary methods provide local argument | 核心特征线与局部Lipschitz/不动点论证移入Methods，删这一句转介 |
| 121 | complete numerical ledger in supporting information/source data | 指向仓库 `research/outputs/refinement.csv`；仍说明有限网格参照不是精确解 |
| 139 | All numerical dose maps in supplement | 完整概率向量指向每个run的 `parameters.json`；正文保留已明确的机会年龄窗与完成恢复定义 |
| 141 | SII multiplication versus unidentified multiplier in supplement | 改为仓库 `research/sources/sii_bridge.csv` 与来源账本；保留原来源未识别限制 |
| 179 | mortality comparison in supplement | 直接指Figure6D及仓库 `research/outputs/mortality_brackets.csv` |
| 210 | all CTMC observation times and interval definitions in supporting information | 指向仓库 `revision/outputs/ctmc_absorbing_event.csv` 与原始replicates；主文已有12月结果与区间类型 |
| 212 | complete gain-refinement ledger in supporting information | 指向仓库 `research/outputs/refinement.csv` |
| 218 | Table aggregation and Figure2 | 如删重复反例表，保留Figure2；案例数、C/H与人数仍在该段和图注，不损失核心证据 |
| 219–228 | aggregation Table1 | 可删，全部数值仍在Figure2/正文与 `research/outputs/aggregation.csv` |
| 233 | trajectories and year10 separately supplied | 明确仓库 `figures/trajectories/` 及 `research/inputs/Final/outputs/stable.csv` |
| 253 | Figure4图注已自洽，无SI引用 | 原样保留即可；完整交互数值由仓库 interactions/shapley/herd CSV承担 |
| 256 | schedule table and Figure5 | 如删重复schedule表，保留Figure5和仓库 `schedules.csv` |
| 257–273 | schedule QALY table | 可删，9种策略×3产品完整数值在 `schedules.csv`，Figure5还呈现病例和成本 |
| 280 | Figure6 and supplementary tables | 指向Figure6与仓库 sensitivity/aom/mortality CSV |
| 301 | migration table isolates accounting changes | 指向仓库 `research/audit/accounting_migration.csv` 和版本说明；主文仍完整定义现行账本 |
| 309 | Supporting information maps all project research, corrections and history | 换成仓库研究记录与结果范围；保留不同估计对象/成本年不混用的科学限定，可删除内部Week名称 |
| 314 | Data/code availability | 明确本次仓库为准备上传的本地release；未发布前不能伪造GitHub URL/DOI/commit号 |
| 321–323 | 注释掉的SI captions | 新稿删除这些失效注释，正式编译包不收SI wrapper、appendix或S1–S9图 |

## 原S1–S15逐节去向

| 原节 | 正文应保留/补入 | 仓库中保留的完整证据 |
|---|---|---|
| S1 特征线与局部适定性 | 特征线公式、正性/平衡条件、比值泛函Lipschitz界和局部结论的作用范围；不把关键证明只移到仓库 | `docs/extended_research_record.md#s1-characteristic-formula-and-local-nonlinear-well-posedness` 全证明记录 |
| S2 有限随机反馈与CTMC分布 | 均值非闭合的简短例子、控制链生成元/初值/时点、终端均值方差/吸收概率定义；现有12月结果可保留 | `research/outputs/ctmc.csv`, `nonlinear_expectation.json`; `revision/outputs/ctmc_variability.csv`, `ctmc_absorbing_event.csv`, `ctmc_absorbing_replicates.csv` |
| S3 数值/经济记账 | 现有正文已含主要数值步骤、incident-flow成本、混合时限和剂次/疾病贴现；只补必要精确定义/变换 | `research/inputs/Final/outputs/runs/*/{parameters.json,economic_parameters.json,economic_details.csv,yearly_unrounded.csv,vaccine_costs.csv}`；全文研究记录 |
| S4 成果覆盖与证据层级 | 正文已说明数学性质/计算验证/参数证据之别，不新增内部成果清单表 | 本映射、`artifact_map.csv`, `docs/data_dictionary.md`, `research/coverage` |
| S5 完整参数表 | 新增关键参数/结构假设表；表中单位和证据状态明确，不称219个常数都已验证 | `research/sources/parameter_sources.csv`, 每run参数JSON，研究记录中的完整两张参数表 |
| S6 闭合队列与AOM分母 | 保留16/16均值检查、AOM端点定义和条件期望型death-reward限制；无须把全部重复模拟表搬入主文 | `research/outputs/cohort/`, `research/inputs/Cohort/` 与历史模型比较数据；研究记录全表 |
| S7 全因子/交互/herd反事实 | 保留Figure4、三产品000/111主结果表、Shapley和herd-off/on定义及关键结论 | `research/inputs/Final/outputs/{factorial,interactions,shapley,herd_contributions}.csv`及37核心run |
| S8 全接种策略 | 保留Figure5、代表性数值、有限机会窗和恢复策略定义；可删重复主表 | `research/inputs/Final/outputs/schedules.csv` 与27 schedule runs |
| S9 全确定性扰动 | 保留Figure6A/B、20设置与PSA区别、成本视角结论 | `research/inputs/Final/outputs/sensitivity.csv` 与对应NV/疫苗runs；完整20设置×3产品 |
| S10 AOM/死亡情景 | 原正文和Figure6C/D已给端点分母、0.303代理、现行死亡分配及其局限，无需额外主表 | `research/outputs/aom.csv`、`mortality_brackets.csv`、`mortality.json`；AOM20行与死亡16行完整保留 |
| S11 全轨迹与年10汇总 | 保留Figure3代表性年龄密度、状态积分、C/r/H、病例轨迹；强调年10不是拟合平衡 | `figures/trajectories/`全9幅；144runs完整states/monthly/herd/balance；`stable.csv`,`trajectories.csv`,`trajectory_diagnostics.csv` |
| S12 全数值账本 | 保留Figure1、最细gain、参照网格差与制造解阶数等现有数值，不新增冗余表 | `research/outputs/`中的refinement、manufactured、ctmc、checks CSV及provenance JSON，与20细化输出；运行验证和源哈希 |
| S13 早期开发验证/弹性/DSA | 主文可用一句说明另有历史队列基准；不把不同成本年、IRR、保护实施或失联PSA汇总混入当前结果 | `research/coverage/inputs`及全部reconstruction outputs；研究记录保留全部数字、证据限定与原始200-draw汇总不可复核事实 |
| S14 参数来源与转换 | 基础年龄风险公式/单位、年龄权重、价格代理与事件相对贴现必要内容进入Methods/关键表；其余完整逐项来源保留仓库 | `research/sources/{parameter_sources,input_comparison,productivity_age_weights,sii_bridge}.csv` 和verification/reference identities；研究记录全文 |
| S15 新旧稿迁移账本 | 保留当前定义和不混用旧分母/成本年的限定；旧稿逐值变化不需主文长表 | `research/audit/{original_core,accounting_migration}.csv`, `research/inputs/Final/outputs/accounting_migration.csv`及研究记录完整历史表 |

## 最小主文资产集合

- 六张主图全部保持现有PDF字节不变。
- 建议主表1为关键输入/假设；主表2为现有三产品000/111病例、死亡reward、QALY、成本、ICER结果，数值和定义不变。
- 原反例表与接种时序表可由已有图/正文和完整CSV承担；删除后需同步改相应两句Table引用。
- 不再有独立Supplement PDF/TeX；仓库研究记录保留原46张显示表及全部限定、102个数学表达式、9个扩展轨迹图入口。
- 不能以“迁至GitHub”替代核心论证，不能把未核实常数改称观测数据，不能把CTMC预测区间或Monte Carlo均值CI改称PSA。
