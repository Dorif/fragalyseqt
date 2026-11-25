# DDEM v0.9 (DNA 数据交换报文) 格式技术规范

## 1. 概述与基本架构

**DDEM (DNA Data Exchange Message) v0.9** 是一种基于 XML 的标准化法医 DNA 数据交换格式，用于在法医物证实验室、国家法医 DNA 数据库管理系统以及实验室信息管理系统 (LIMS) 之间传输、共享和批量导入 STR (短串联重复序列) 基因分型数据。

从底层数据架构上看，DDEM v0.9 采用并裁剪了国家 DNA 数据库导入格式规范（命名空间 `urn:NDNADImportFile-schema`），并针对法医物证检案规范进行了严格的本地化与结构约束：
- **XML 命名空间**：`urn:NDNADImportFile-schema`
- **模式校验标准**：W3C XML Schema 1.0
- **根元素**：`<NDNADImportFile>`
- **字符编码**：UTF-8
- **校验原则**：严格命名空间限定 (`elementFormDefault="qualified"`)、限定长度的简单类型定义、枚举约束以及文档级唯一性约束。

---

## 2. 格式结构比对：DDEM v0.9 与 CODIS CMF 3.2

尽管两者均用于封装法医 STR 分型图谱数据，但 DDEM v0.9 在层级结构、元数据挂载点以及术语定义上与通用 CODIS CMF 3.2 (Combined DNA Index System Common Message Format) 存在显著差异：

| 特性 / 字段 | CODIS CMF 3.2 | DDEM v0.9 | 结构与架构差异解析 |
| :--- | :--- | :--- | :--- |
| **根节点** | `<CODISImportFile>` | `<NDNADImportFile>` | 根容器采用国家导入格式标准命名空间。 |
| **目标命名空间** | CODIS 专用命名空间或缺省 | `urn:NDNADImportFile-schema` | 全文元素必须由目标命名空间完全限定。 |
| **文件头版本** | `<HEADERVERSION>3.2</HEADERVERSION>` | `<HEADERVERSION>0.9</HEADERVERSION>` | 数值类型为 `xsd:decimal`（总长度 2 位，小数位 1 位）。 |
| **报文类型** | `<MESSAGETYPE>` | `<MESSAGETYPE>Import</MESSAGETYPE>` | 严格限定为字符串枚举值 `"Import"`（固定 6 字符）。 |
| **接收实验室** | `<DESTINATIONORI>` | `<DESTINATIONLAB>` | 接收方实验室代码标签（最大长度 10 字符）。 |
| **报送实验室** | `<SOURCEORI>` | `<SOURCELAB>` | 发送方实验室代码标签（最大长度 10 字符）。 |
| **样本类别分类** | 英文自由文本 / CODIS 标准类别 | 严格限定 22 种规范汉字类别 (`SpecimenCategoryType`) | 必须使用标准简体中文（如 `嫌疑人`、`受害人`、`未知样本` 等）。 |
| **试剂盒与批次信息** | 通常位于文件级或样本级 | 在文件头支持，并可在基因座 (`LOCUS`) 属性中指定 | 属性 `KIT` 与 `BATCHID` 直接下沉绑定至各个基因座级别。 |
| **等位基因匹配属性** | 不强制或使用特定标记 | 显式属性 `ALLELEREQUIRED="false"` / `"true"` | 每个 `<ALLELE>` 标签带有布尔匹配控制属性。 |
| **基因座唯一性约束** | 业务逻辑层面控制 | XML Schema 强制：`<xsd:unique name="UNIQUE_LOCI">` | 单个样本内严禁出现重复基因座。 |
| **样本编号唯一性约束** | 数据库系统层面控制 | XML Schema 强制：`<xsd:unique name="UNIQUE_SPEC">` | 单个 XML 文件内严禁出现重复的 `<SPECIMENID>`。 |

---

## 3. 核心 XML 元素层级与数据字典规范

根据 DDEM v0.9 模式定义，文档层级分为三层：报文头信息、样本信息 (`SPECIMEN`)、基因座与等位基因信息 (`LOCUS` / `ALLELE`)。

### 3.1. 报文头元素 (`<NDNADImportFile>`)

根节点包含了文件报送批次信息及样本集：

```xml
<NDNADImportFile xmlns="urn:NDNADImportFile-schema">
    <HEADERVERSION>0.9</HEADERVERSION>
    <MESSAGETYPE>Import</MESSAGETYPE>
    <DESTINATIONLAB>LAB_NAT01</DESTINATIONLAB>
    <SOURCELAB>LAB_PROV08</SOURCELAB>
    <SUBMITBYUSERID>analyst_wang</SUBMITBYUSERID>
    <SUBMITDATETIME>2026-11-25T14:30:00</SUBMITDATETIME>
    <BATCHID>BATCH_20261125_01</BATCHID>
    <KIT>Identifiler</KIT>
    <!-- 一个或多个 SPECIMEN 元素 -->
</NDNADImportFile>
```

#### 报文头元素清单：
1. **`<HEADERVERSION>`** (*可选, minOccurs=0, maxOccurs=1*):
   - 类型：`DDEMHeaderVersionType` (`xsd:decimal`，总位数 2，小数位数 1)。
   - 当前固定取值：`0.9`。
2. **`<MESSAGETYPE>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`DDEMMessageType` (`xsd:string`，固定长度 6)。
   - 唯一有效枚举取值：`Import`。
3. **`<DESTINATIONLAB>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`DDEMLabType` (`xsd:string`，长度 1～10 字符)。
   - 接收数据包的目标实验室或数据库节点编号。
4. **`<SOURCELAB>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`DDEMLabType` (`xsd:string`，长度 1～10 字符)。
   - 报送数据的源实验室编号。
5. **`<SUBMITBYUSERID>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`DDEMUserType` (`xsd:string`，长度 1～20 字符)。
   - 提交报文的检验人员或操作员登录账号。
6. **`<SUBMITDATETIME>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`NDNADImportDate` (`xsd:dateTime`)。
   - 格式遵循 ISO-8601 (`YYYY-MM-DDThh:mm:ss`)。有效区间：`1900-01-01T00:00:00` 至 `2079-06-06T00:00:00`。
7. **`<BATCHID>`** (*可选, minOccurs=0, maxOccurs=1*):
   - 类型：`BatchIDType` (`xsd:string`，长度 0～32 字符)。
   - 扩增或检测批次/孔板流水标识。
8. **`<KIT>`** (*可选, minOccurs=0, maxOccurs=1*):
   - 类型：`KitType` (`xsd:string`，长度 0～32 字符)。
   - 主 STR 荧光检测试剂盒名称（必须匹配预定义枚举）。
### 3.2. 样本元素 (`<SPECIMEN>`)

每个 `<SPECIMEN>` 元素（minOccurs=1, maxOccurs=unbounded）封装了一个独立的法医 STR DNA 基因分型图谱。

```xml
<SPECIMEN CASEID="CASE-2026-0812" SOURCEID="Yes" PARTIAL="false">
    <SPECIMENID>SPEC-2026-0091</SPECIMENID>
    <SPECIMENCATEGORY>嫌疑人</SPECIMENCATEGORY>
    <SPECIMENCOMMENT>现场擦拭物检材</SPECIMENCOMMENT>
    <!-- 1 到 40 个 LOCUS 元素 -->
</SPECIMEN>
```

#### `<SPECIMEN>` 属性定义：
- **`CASEID`** (*可选, `CaseIDType`*): 案件编号字符串（长度 0～32 字符）。
- **`SOURCEID`** (*可选, `SourceIDType`*): 枚举字符串（取值为 `"Yes"`、`"No"`、`"N/A"`），标明检材来源个体身份是否明确。
- **`PARTIAL`** (*可选, `xsd:boolean`*): 布尔值（`true` / `false` / `1` / `0`），标明当前图谱是否为部分分型（降解检材）。

#### `<SPECIMEN>` 子元素清单：
1. **`<SPECIMENID>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`SpecimenIDType` (`xsd:string`，长度 1～24 字符)。
   - 样本唯一识别码。XML 模式通过 `<xsd:unique name="UNIQUE_SPEC">` 强制保证单份文件内该标识唯一。
2. **`<SPECIMENCATEGORY>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`SpecimenCategoryType` (`xsd:string`，长度 1～21 字符)。
   - 样本分类。必须严格取自国家规范规定的 22 种规范汉字类别之一（见第 4.1 节）。
3. **`<SPECIMENCOMMENT>`** (*可选, minOccurs=0, maxOccurs=1*):
   - 类型：`SpecimenCommentType` (`xsd:string`，长度 0～255 字符)。
   - 检验人员填写的样本备注、检材来源或实验说明。
4. **`<LOCUS>`** (*必填, minOccurs=1, maxOccurs=40*):
   - 基因座分型信息。

---

### 3.3. 基因座元素 (`<LOCUS>`)

每个样本下包含 1 到 40 个 `<LOCUS>` 元素：

```xml
<LOCUS BATCHID="BATCH_20261125_01" KIT="Identifiler">
    <LOCUSNAME>D8S1179</LOCUSNAME>
    <READINGBY>analyst_wang</READINGBY>
    <READINGDATETIME>2026-11-25T14:35:00</READINGDATETIME>
    <ALLELE ALLELEREQUIRED="false">
        <ALLELEVALUE>13</ALLELEVALUE>
    </ALLELE>
    <ALLELE ALLELEREQUIRED="false">
        <ALLELEVALUE>15</ALLELEVALUE>
    </ALLELE>
</LOCUS>
```

#### `<LOCUS>` 属性定义：
- **`BATCHID`** (*可选, `BatchIDType`*): 扩增检测批次/PCR 孔板编号（长度 0～32 字符）。
- **`KIT`** (*可选, `KitType`*): 该基因座所属的荧光复合扩增检测试剂盒名称。

#### `<LOCUS>` 子元素清单：
1. **`<LOCUSNAME>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`LocusNameType` (`xsd:string`，长度 1～10 字符)。
   - 基因座名称，必须严格属于 22 个受支持标记枚举之一（见第 4.3 节）。
   - 模式通过 `<xsd:unique name="UNIQUE_LOCI">` 约束单个样本内基因座名称严禁重复。
2. **`<READINGBY>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`DDEMUserType` (`xsd:string`，长度 1～20 字符)。
   - 审核或确认该基因座分型结果的法医专业技术人员账号。
3. **`<READINGDATETIME>`** (*必填, minOccurs=1, maxOccurs=1*):
   - 类型：`NDNADImportDate` (`xsd:dateTime`)。
   - 判读确认时间戳（ISO-8601 格式：`YYYY-MM-DDThh:mm:ss`）。
4. **`<ALLELE>`** (*必填, minOccurs=1, maxOccurs=8*):
   - 等位基因分型数值容器（支持 1～8 个等位基因，兼容单倍体、纯合子、杂合子、复等位基因及复杂混合分型）。

---

### 3.4. 等位基因元素 (`<ALLELE>`)

```xml
<ALLELE ALLELEREQUIRED="false">
    <ALLELEVALUE>12</ALLELEVALUE>
</ALLELE>
```

- **属性 `ALLELEREQUIRED`** (*可选, 缺省值为 `false`, `xsd:boolean`*):
  - 布尔标志（`true` / `false` / `1` / `0`）。
  - 指示在数据库比对与串并案检索中该等位基因是否为强制必配等位基因（常用于降解检材筛选或混合物拆分）。
- **子元素 `<ALLELEVALUE>`** (*必填, minOccurs=1, maxOccurs=1*):
  - 类型：`AlleleValueType` (`xsd:string`，长度 1～10 字符)。
  - 等位基因呼叫值（例如：`9.3`、`14`、`X`、`Y`、`OL`、`INC` 等）。

---

## 4. 枚举字典与受限数据类型规范

DDEM v0.9 采用封闭式白名单设计，所有样本分类、检测试剂盒及基因座名称必须与预定义枚举完全一致。

### 4.1. 样本分类枚举 (`SpecimenCategoryType`)

模式定义了 22 种规范简体汉字样本类别（最大长度 21 字符）：

| 序号 | 规范中文枚举名称 (XML 取值) | 拼音 | 类别说明 | 对应国际通用分类概念 |
| :---: | :--- | :--- | :--- | :--- |
| 1 | `嫌疑父` | Xiányí fù | 亲权鉴定中被指认的争议父亲 | Alleged Father |
| 2 | `嫌疑母` | Xiányí mǔ | 亲权鉴定中被指认的争议母亲 | Alleged Mother |
| 3 | `子女` | Zǐnǚ | 被鉴定的亲生子女（子或女） | Biological Child |
| 4 | `生父` | Shēngfù | 经鉴定的确认亲生父亲 | Biological Father |
| 5 | `生母` | Shēngmǔ | 经鉴定的确认亲生母亲 | Biological Mother |
| 6 | `同胞` | Tóngbāo | 亲兄弟姐妹（全同胞或半同胞） | Biological Sibling |
| 7 | `群体遗传学样本` | Qúntǐ yíchuánxué yàngběn | 人群基因多态性调查与频率统计样本 | Population Study Sample |
| 8 | `能力测试样本` | Nénglì cèshì yàngběn | 实验室能力验证、盲样考核样本 | Proficiency Test Sample |
| 9 | `父系样本` | Fùxì yàngběn | 具有明确父系血缘的亲属样本 | Paternal Relative |
| 10 | `未知样本` | Wèizhī yàngběn | 案件现场物证提取的未知源 STR 检材 | Forensic Unknown / Crime Scene Evidence |
| 11 | `排除嫌疑个体` | Páichú xiányí gètǐ | 经核查排除涉案嫌疑的对照个体 | Elimination Reference |
| 12 | `罪犯` | Zuìfàn | 依法入库的已决犯与服刑人员 | Convicted Offender |
| 13 | `失踪人员` | Shīzōng rényuán | 失踪或被拐卖人员生前物证或参考样本 | Missing Person |
| 14 | `混合样本` | Hùnhé yàngběn | 包含两个或以上个体成分的混合 DNA 检材 | Forensic Mixture |
| 15 | `母系样本` | Mǔxì yàngběn | 具有明确母系血缘的亲属样本 | Maternal Relative |
| 16 | `双亲关系` | Shuāngqīn guānxì | 包含双亲信息的亲子鉴定参照关系 | Biological Parents (Duo/Trio) |
| 17 | `受害人` | Shòuhàirén | 刑事案件被害人已知人身比对样本 | Victim, Known |
| 18 | `已死亡个体` | Yǐ sǐwáng gètǐ | 无名尸体、遇难者遗骸或骨骼检材 | Deceased Individual / Unidentified Remains |
| 19 | `嫌疑人` | Xiányí rén | 涉案刑事侦查在押或排查嫌疑人员 | Suspect, Known |
| 20 | `质控样本` | Zhìkòng yàngběn | 阳性对照 (9947A/2800M)、阴性对照或等位基因阶梯 | Quality Control / Ladder |
| 21 | `配偶` | Pèi'ǒu | 婚姻配偶对照样本 | Spouse |
| 22 | `其他` | Qítā | 无法归入上述类别的其他法医样本 | Other |

---

### 4.2. 复合扩增检测试剂盒枚举 (`KitType`)

模式定义了 15 个受支持的 PCR 试剂盒标识（最大长度 32 字符）。部分枚举值包含固有的前导或后随空格，生成 XML 时必须严格保留该字符串格式：

| 序号 | 试剂盒枚举字面值 (严格匹配) | 试剂盒名称与特征 | 生产商 / 研发背景 |
| :---: | :--- | :--- | :--- |
| 1 | `GlobalFiler` | GlobalFiler 24 基因座荧光检测系统 | Thermo Fisher Scientific |
| 2 | `Identifiler` | AmpFLSTR Identifiler 16 基因座系统 | Thermo Fisher Scientific |
| 3 | `Profiler Plus` | AmpFLSTR Profiler Plus 系统 | Thermo Fisher Scientific |
| 4 | `PowerPlex 18` | PowerPlex 18D / 18 系统 | Promega Corporation |
| 5 | `PowerPlex 21` | PowerPlex 21 系统 (21 位点) | Promega Corporation |
| 6 | `PowerPlex 16` | PowerPlex 16 系统 | Promega Corporation |
| 7 | `DNATyper 15` | DNATyper 15 中国法医标准 15 位点系统 | 公安部物证鉴定中心 (MPS PRC) |
| 8 | `Sinofiler` | AmpFLSTR Sinofiler 系统 (针对中国群体优化) | Applied Biosystems |
| 9 | `SGM Plus` | AmpFLSTR SGM Plus 系统 | Applied Biosystems |
| 10 | ` ACGU 17+1` *(注意：带 1 个前导空格)* | ACGU 17+1 荧光检测试剂盒 | ACGU (中国) |
| 11 | ` GoldenEye 16A ` *(注意：前后各 1 个空格)* | GoldenEye DNA 身份鉴定系统 16A | Peoplespot 基点认知 (中国) |
| 12 | ` DNATyper 15 Plus ` *(注意：前后各 1 个空格)* | DNATyper 15 Plus 系统 | 公安部物证鉴定中心 (MPS PRC) |
| 13 | ` GoldenEye 16BT` *(注意：带 1 个前导空格)* | GoldenEye 16BT 系统 | Peoplespot 基点认知 (中国) |
| 14 | ` GoldenEye 16C` *(注意：带 1 个前导空格)* | GoldenEye 16C 系统 | Peoplespot 基点认知 (中国) |
| 15 | ` GoldenEye 20A ` *(注意：前后各 1 个空格)* | GoldenEye 20A 系统 (20 位点) | Peoplespot 基点认知 (中国) |

---

### 4.3. 基因座名称枚举 (`LocusNameType`)

DDEM v0.9 严格收录并限定了 22 个常染色体 STR 及性别鉴定标记：

```text
AMEL, Amelogenin, CSF1PO, D13S317, D16S539, D18S51, D19S433, D21S11,
D2S1338, D3S1358, D5S818, D7S820, D8S1179, FGA, Penta D, Penta E,
TH01, TPOX, vWA, D22S1045, D2S441
```

*说明：性别决定位点在模式中支持两种写法：`AMEL` 与 `Amelogenin`。单个样本中仅允许出现其中一种。*

---

## 5. 数据完整性与模式校验规则

生成符合 DDEM v0.9 规范的 XML 报文必须满足以下 W3C Schema 约束：

1. **目标命名空间一致性**：根元素必须声明 `xmlns="urn:NDNADImportFile-schema"`，文档中所有子元素均属于该命名空间。
2. **文档级唯一性约束**：
   - `<xsd:unique name="UNIQUE_SPEC">`：单份 XML 报文内，所有样本的 `<SPECIMENID>` 必须唯一。
   - `<xsd:unique name="UNIQUE_LOCI">`：单个 `<SPECIMEN>` 内，所有基因座的 `<LOCUSNAME>` 必须唯一。
3. **数量上限与下限**：
   - `<SPECIMEN>`：至少包含 1 个样本，上限不限 (`minOccurs="1"`, `maxOccurs="unbounded"`).
   - `<LOCUS>`：单个样本内包含 1 到 40 个基因座 (`minOccurs="1"`, `maxOccurs="40"`).
   - `<ALLELE>`：单个基因座内包含 1 到 8 个等位基因 (`minOccurs="1"`, `maxOccurs="8"`).
4. **字符串长度与时间区间约束**：
   - `DESTINATIONLAB` / `SOURCELAB`: 1～10 字符。
   - `SPECIMENID`: 1～24 字符。
   - `SUBMITBYUSERID` / `READINGBY`: 1～20 字符。
   - `ALLELEVALUE`: 1～10 字符。
   - `CASEID` / `BATCHID`: 0～32 字符。
   - `SPECIMENCOMMENT`: 0～255 字符。
   - 时间戳范围：`1900-01-01T00:00:00` 至 `2079-06-06T00:00:00`。
5. **元素顺序严格性**：文档内子元素出现顺序必须严格按照 XSD `<xsd:sequence>` 定义的排列顺序输出。

---

## 6. DDEM v0.9 标准 XML 实例

```xml
<?xml version="1.0" encoding="UTF-8"?>
<NDNADImportFile xmlns="urn:NDNADImportFile-schema">
    <HEADERVERSION>0.9</HEADERVERSION>
    <MESSAGETYPE>Import</MESSAGETYPE>
    <DESTINATIONLAB>LAB_NAT01</DESTINATIONLAB>
    <SOURCELAB>LAB_PROV08</SOURCELAB>
    <SUBMITBYUSERID>operator_01</SUBMITBYUSERID>
    <SUBMITDATETIME>2026-11-25T14:30:00</SUBMITDATETIME>
    <BATCHID>PLATE_20261125_STR</BATCHID>
    <KIT>Identifiler</KIT>
    <SPECIMEN CASEID="CASE_2026_0981" SOURCEID="Yes" PARTIAL="false">
        <SPECIMENID>SAMPLE_88291</SPECIMENID>
        <SPECIMENCATEGORY>嫌疑人</SPECIMENCATEGORY>
        <SPECIMENCOMMENT>现场擦拭物检材</SPECIMENCOMMENT>
        <LOCUS BATCHID="PLATE_20261125_STR" KIT="Identifiler">
            <LOCUSNAME>Amelogenin</LOCUSNAME>
            <READINGBY>operator_01</READINGBY>
            <READINGDATETIME>2026-11-25T15:00:00</READINGDATETIME>
            <ALLELE ALLELEREQUIRED="false">
                <ALLELEVALUE>X</ALLELEVALUE>
            </ALLELE>
            <ALLELE ALLELEREQUIRED="false">
                <ALLELEVALUE>Y</ALLELEVALUE>
            </ALLELE>
        </LOCUS>
        <LOCUS BATCHID="PLATE_20261125_STR" KIT="Identifiler">
            <LOCUSNAME>D8S1179</LOCUSNAME>
            <READINGBY>operator_01</READINGBY>
            <READINGDATETIME>2026-11-25T15:00:00</READINGDATETIME>
            <ALLELE ALLELEREQUIRED="false">
                <ALLELEVALUE>12</ALLELEVALUE>
            </ALLELE>
            <ALLELE ALLELEREQUIRED="false">
                <ALLELEVALUE>14</ALLELEVALUE>
            </ALLELE>
        </LOCUS>
    </SPECIMEN>
</NDNADImportFile>
```
