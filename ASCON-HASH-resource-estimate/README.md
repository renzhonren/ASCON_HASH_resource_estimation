# ASCON-HASH ProjectQ 资源估算

主程序逐门构建论文中的 QAND 异地、Amy 异地浅层和 Amy 原地三个方案，输出量子比特数、总深度及基本门数。

## 准备线性扩散门序

主程序使用 [Roy 等人的 ASCON 线性层门序](https://github.com/sohamroy19/ascon-linear-layer/blob/main/implementations/xzlbz.txt)，但本目录不随附该文件。运行前请从原作者仓库下载，在本目录保存为 `roy_linear_xzlbz.txt`。PowerShell 命令：

```powershell
Invoke-WebRequest -Uri 'https://github.com/sohamroy19/ascon-linear-layer/raw/refs/heads/main/implementations/xzlbz.txt' -OutFile 'roy_linear_xzlbz.txt'
```

## 安装与运行

在本目录执行：

```powershell
python -m pip install -r requirements.txt
python ASCON_HASH_resource_estimation.py --scheme all
```

`--scheme` 可选 `qand`、`amy-outplace`、`amy-inplace` 或 `all`。

`qubits` 字段排除 256 比特消息输入和 256 比特摘要寄存器；`physical_peak_qubits_including_message_and_digest` 字段给出包含两者的实际峰值。

`--include-fixed-initial-state` 可额外计入固定初态制备；异地方案只执行后续复用所需的逆轮：QAND 异地为 78 次，Amy 异地为 77 次。程序使用 ProjectQ `ResourceCounter.depth_of_dag` 计算完整门序的总深度。

| 默认方案 | 论文口径量子比特 | 实际峰值量子比特 | 总深度 | T 门 | CNOT 门 | 单比特 Clifford 门 | 测量 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| QAND 异地 | 6,720 | 7,232 | 21,233 | 122,880 | 784,666 | 204,194 | 24,960 |
| Amy 异地浅层 | 7,040 | 7,552 | 12,681 | 387,520 | 843,167 | 123,166 | 0 |
| Amy 原地 | 320 | 832 | 17,097 | 301,056 | 448,544 | 166,273 | 0 |


