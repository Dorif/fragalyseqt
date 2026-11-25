# DDEM v0.9 (DNA Data Exchange Message) Format Specification

## 1. Overview and Architecture

**DDEM (DNA Data Exchange Message) v0.9** is an XML-based data exchange format designed for transferring and importing forensic Short Tandem Repeat (STR) DNA profiles into national forensic DNA database management systems and LIMS platforms.

Architecturally, DDEM v0.9 is structured around the United Kingdom National DNA Database (UK NDNAD) Import File Schema (`urn:NDNADImportFile-schema`), adapted and constrained for Chinese forensic casework standards and national DNA database requirements:
- **XML Namespace**: `urn:NDNADImportFile-schema`
- **Validation Schema**: W3C XML Schema 1.0
- **Root Element**: `<NDNADImportFile>`
- **Character Encoding**: UTF-8
- **Validation**: Strict schema validation with qualified element names (`elementFormDefault="qualified"`), explicit simple types, field length restrictions, and document-level uniqueness constraints.

---

## 2. Structural Differences: DDEM v0.9 vs. CODIS CMF 3.2

While both formats encapsulate forensic STR typing profiles, DDEM v0.9 differs significantly in document hierarchy, metadata placement, and taxonomy compared to standard CODIS CMF 3.2 (Combined DNA Index System Common Message Format):

| Feature / Element | CODIS CMF 3.2 | DDEM v0.9 | Architectural Difference |
| :--- | :--- | :--- | :--- |
| **Root Element** | `<CODISImportFile>` | `<NDNADImportFile>` | Root container aligned with the NDNAD schema specification. |
| **Target Namespace** | Default or CODIS schema namespace | `urn:NDNADImportFile-schema` | Explicit namespace qualification required on all elements. |
| **Header Version** | `<HEADERVERSION>3.2</HEADERVERSION>` | `<HEADERVERSION>0.9</HEADERVERSION>` | Numeric format constrained to `xsd:decimal` (totalDigits=2, fractionDigits=1). |
| **Message Type** | `<MESSAGETYPE>` | `<MESSAGETYPE>Import</MESSAGETYPE>` | Strictly constrained to enumeration `"Import"`. |
| **Destination Lab Tag** | `<DESTINATIONORI>` | `<DESTINATIONLAB>` | Renamed element representing the recipient laboratory identifier (max 10 chars). |
| **Source Lab Tag** | `<SOURCEORI>` | `<SOURCELAB>` | Renamed element representing the sender laboratory identifier (max 10 chars). |
| **Specimen Category** | Free-text / Standard CODIS Categories (e.g., `Forensic Unknown`, `Convicted Offender`) | Strict enumeration of 22 Chinese forensic categories (`SpecimenCategoryType`) | Requires exact Chinese Simplified characters (e.g., `嫌疑人`, `受害人`, `未知样本`). |
| **Batch / Kit Metadata** | Often placed at file or specimen level | Supported at header level and explicitly assigned per locus via XML attributes | `<LOCUS BATCHID="..." KIT="...">` embeds run metadata directly into each locus element. |
| **Locus Name List** | Broad list of international STR loci | Strict enumeration of 22 supported loci (`LocusNameType`) | Limited to predefined national database markers. |
| **PCR Kit List** | Flexible string | Strict enumeration of 15 commercial STR multiplex kits (`KitType`) | Constrained to validated domestic and international STR systems. |
| **Allele Required Flag** | Not present | Attribute `ALLELEREQUIRED` on `<ALLELE>` (default: `false`) | Boolean flag indicating whether the allele is obligate during database searching. |
| **Uniqueness Constraints** | Handled at database ingestion level | Enforced by `<xsd:unique>` in XML schema | Schema guarantees unique `SPECIMENID` per file and unique `LOCUSNAME` per specimen. |

---

## 3. Detailed Data Dictionary and Hierarchy

```
<NDNADImportFile xmlns="urn:NDNADImportFile-schema">
  ├── <HEADERVERSION>0.9</HEADERVERSION>
  ├── <MESSAGETYPE>Import</MESSAGETYPE>
  ├── <DESTINATIONLAB>...</DESTINATIONLAB>
  ├── <SOURCELAB>...</SOURCELAB>
  ├── <SUBMITBYUSERID>...</SUBMITBYUSERID>
  ├── <SUBMITDATETIME>...</SUBMITDATETIME>
  ├── <BATCHID>...</BATCHID>                     [optional]
  ├── <KIT>...</KIT>                             [optional]
  └── <SPECIMEN ...>                             [1..*]
        ├── <SPECIMENID>...</SPECIMENID>
        ├── <SPECIMENCATEGORY>...</SPECIMENCATEGORY>
        ├── <SPECIMENCOMMENT>...</SPECIMENCOMMENT> [optional]
        └── <LOCUS ...>                          [1..40]
              ├── <LOCUSNAME>...</LOCUSNAME>
              ├── <READINGBY>...</READINGBY>
              ├── <READINGDATETIME>...</READINGDATETIME>
              └── <ALLELE ...>                   [1..8]
                    └── <ALLELEVALUE>...</ALLELEVALUE>
```

### 3.1. File Header Elements (`<NDNADImportFile>`)

The document root `<NDNADImportFile>` mandates an exact sequence of header tags followed by one or more specimen elements:

1. **`<HEADERVERSION>`** *(Optional, minOccurs=0, maxOccurs=1)*:
   - Type: `DDEMHeaderVersionType` (`xsd:decimal`, max 2 digits, 1 decimal place).
   - Expected value: `0.9`.
2. **`<MESSAGETYPE>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `DDEMMessageType` (`xsd:string`, fixed length 6).
   - Permitted enumeration: `"Import"`.
3. **`<DESTINATIONLAB>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `DDEMLabType` (`xsd:string`, length 1 to 10 chars).
   - Alphanumeric code identifying the receiving laboratory or database node.
4. **`<SOURCELAB>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `DDEMLabType` (`xsd:string`, length 1 to 10 chars).
   - Alphanumeric code identifying the submitting laboratory or agency.
5. **`<SUBMITBYUSERID>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `DDEMUserType` (`xsd:string`, length 1 to 20 chars).
   - User account or identifier of the analyst/operator submitting the batch.
6. **`<SUBMITDATETIME>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `NDNADImportDate` (`xsd:dateTime`, ISO-8601 format: `YYYY-MM-DDThh:mm:ss`).
   - Must fall within the range `1900-01-01T00:00:00` to `2079-06-06T00:00:00`.
7. **`<BATCHID>`** *(Optional, minOccurs=0, maxOccurs=1)*:
   - Type: `BatchIDType` (`xsd:string`, length 0 to 32 chars).
   - Overall analytical batch, run plate, or processing session identifier.
8. **`<KIT>`** *(Optional, minOccurs=0, maxOccurs=1)*:
   - Type: `KitType` (`xsd:string`, length 0 to 32 chars).
   - Default STR amplification kit identifier, selected from the enumerated `KitType` list.
### 3.2. Specimen Element (`<SPECIMEN>`)

Each `<SPECIMEN>` element (minOccurs=1, maxOccurs=unbounded) encapsulates an individual STR DNA profile.

```xml
<SPECIMEN CASEID="CASE-2026-0812" SOURCEID="Yes" PARTIAL="false">
    <SPECIMENID>SPEC-2026-0091</SPECIMENID>
    <SPECIMENCATEGORY>嫌疑人</SPECIMENCATEGORY>
    <SPECIMENCOMMENT>Forensic evidence swab</SPECIMENCOMMENT>
    <!-- 1 to 40 LOCUS elements -->
</SPECIMEN>
```

#### Attributes of `<SPECIMEN>`:
- **`CASEID`** *(Optional, `CaseIDType`)*: Alphanumeric case identifier string (max 32 characters).
- **`SOURCEID`** *(Optional, `SourceIDType`)*: Enumeration (`"Yes"`, `"No"`, `"N/A"`), indicating whether the donor source is known.
- **`PARTIAL`** *(Optional, `xsd:boolean`)*: Boolean flag (`true` / `false` / `1` / `0`) indicating whether the genetic profile is partial.

#### Child Elements of `<SPECIMEN>`:
1. **`<SPECIMENID>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `SpecimenIDType` (`xsd:string`, length 1 to 24 chars).
   - Uniquely identifies the specimen within the file (enforced by `<xsd:unique name="UNIQUE_SPEC">`).
2. **`<SPECIMENCATEGORY>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `SpecimenCategoryType` (`xsd:string`, length 1 to 21 chars).
   - Restricted to exactly 22 standardized Chinese forensic categories (see Section 4.1).
3. **`<SPECIMENCOMMENT>`** *(Optional, minOccurs=0, maxOccurs=1)*:
   - Type: `SpecimenCommentType` (`xsd:string`, length 0 to 255 chars).
   - Free-text commentary, sample description, or bench notes.
4. **`<LOCUS>`** *(Mandatory, minOccurs=1, maxOccurs=40)*:
   - Contains locus typing and allele call information.

---

### 3.3. Locus Element (`<LOCUS>`)

Each specimen contains between 1 and 40 `<LOCUS>` elements:

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

#### Attributes of `<LOCUS>`:
- **`BATCHID`** *(Optional, `BatchIDType`)*: Specific processing batch or analytical plate ID (max 32 chars).
- **`KIT`** *(Optional, `KitType`)*: STR amplification kit used for this specific locus (selected from `KitType` enumeration).

#### Child Elements of `<LOCUS>`:
1. **`<LOCUSNAME>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `LocusNameType` (`xsd:string`, length 1 to 10 chars).
   - Must match one of the recognized loci defined in the schema (see Section 4.2). Locus names cannot repeat within a single specimen (enforced by `<xsd:unique name="UNIQUE_LOCI">`).
2. **`<READINGBY>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `DDEMUserType` (`xsd:string`, length 1 to 20 chars).
   - Identifier of the analyst who confirmed or reviewed this locus call.
3. **`<READINGDATETIME>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
   - Type: `NDNADImportDate` (`xsd:dateTime`, `YYYY-MM-DDThh:mm:ss`).
   - Timestamp when the locus interpretation was finalized.
4. **`<ALLELE>`** *(Mandatory, minOccurs=1, maxOccurs=8)*:
   - Sequence of between 1 and 8 allele calls for the locus.

---

### 3.4. Allele Element (`<ALLELE>`)

```xml
<ALLELE ALLELEREQUIRED="false">
    <ALLELEVALUE>12</ALLELEVALUE>
</ALLELE>
```

- **Attribute `ALLELEREQUIRED`** *(Optional, `xsd:boolean`, default: `false`)*:
  - Specifies whether this allele is obligatory during database searching and matching algorithms.
- **Child Element `<ALLELEVALUE>`** *(Mandatory, minOccurs=1, maxOccurs=1)*:
  - Type: `AlleleValueType` (`xsd:string`, length 1 to 10 chars).
  - Alphanumeric allele call designation (e.g., `9.3`, `14`, `X`, `Y`, `OL`, `INC`).

---

## 4. Standard Data Dictionary and Taxonomies

### 4.1. Specimen Categories (`SpecimenCategoryType`)

DDEM v0.9 strictly restricts `<SPECIMENCATEGORY>` to 22 standardized Simplified Chinese terms:

| # | Chinese Term (PRC) | Pinyin | English Translation | Equivalent CODIS Category |
|---|:---|:---|:---|:---|
| 1 | `嫌疑父` | Xiányí fù | Alleged Father | Alleged Father |
| 2 | `嫌疑母` | Xiányí mǔ | Alleged Mother | Alleged Mother |
| 3 | `子女` | Zǐnǚ | Biological Child | Biological Child |
| 4 | `生父` | Shēngfù | Biological Father | Biological Father |
| 5 | `生母` | Shēngmǔ | Biological Mother | Biological Mother |
| 6 | `同胞` | Tóngbāo | Biological Sibling | Biological Sibling |
| 7 | `群体遗传学样本` | Qúntǐ yíchuánxué yàngběn | Population Genetics Sample | Population |
| 8 | `能力测试样本` | Nénglì cèshì yàngběn | Proficiency Testing Sample | Proficiency |
| 9 | `父系样本` | Fùxì yàngběn | Paternal Lineage Sample | Paternal Relative |
| 10 | `未知样本` | Wèizhī yàngběn | Forensic Unknown / Crime Scene Profile | Forensic Unknown |
| 11 | `排除嫌疑个体` | Páichú xiányí gètǐ | Elimination Individual | Elimination, Known |
| 12 | `罪犯` | Zuìfàn | Convicted Offender | Convicted Offender |
| 13 | `失踪人员` | Shīzōng rényuán | Missing Person | Missing Person |
| 14 | `混合样本` | Hùnhé yàngběn | Forensic Mixture | Forensic Mixture |
| 15 | `母系样本` | Mǔxì yàngběn | Maternal Lineage Sample | Maternal Relative |
| 16 | `双亲关系` | Shuāngqīn guānxì | Biological Parents (Trio/Duet) | Biological Parent |
| 17 | `受害人` | Shòuhàirén | Victim | Victim, Known |
| 18 | `已死亡个体` | Yǐ sǐwáng gètǐ | Deceased Individual / Remains | Unidentified Person / Remains |
| 19 | `嫌疑人` | Xiányí rén | Suspect | Suspect, Known |
| 20 | `质控样本` | Zhìkòng yàngběn | Quality Control Sample | Quality Control |
| 21 | `配偶` | Pèi'ǒu | Spouse | Spouse |
| 22 | `其他` | Qítā | Other | Other |

---

### 4.2. Supported STR Loci (`LocusNameType`)

The DDEM v0.9 schema specification defines 22 valid enumeration strings for `LocusNameType`:

1. `AMEL`
2. `Amelogenin`
3. `CSF1PO`
4. `D13S317`
5. `D16S539`
6. `D18S51`
7. `D19S433`
8. `D21S11`
9. `D2S1338`
10. `D3S1358`
11. `D5S818`
12. `D7S820`
13. `D8S1179`
14. `FGA`
15. `Penta D`
16. `Penta E`
17. `TH01`
18. `TPOX`
19. `vWA`
20. `D22S1045`
21. `D2S441`

*(Note: While some companion reference tables list `D6S1043`, the DDEM v0.9 schema strict enumeration restricts the valid set to the 22 names above. Inclusion of unlisted markers will fail schema validation).*

---

### 4.3. Supported PCR Amplification Kits (`KitType`)

The `KitType` simple type defines exact string enumerations for commercial STR typing kits. **String padding (leading/trailing whitespace) present in the XSD specification must be preserved verbatim**:

| # | Schema Enumeration Value (`KitType`) | Kit Manufacturer / Origin |
|---|:---|:---|
| 1 | `"GlobalFiler"` | Thermo Fisher Scientific |
| 2 | `"Identifiler"` | Thermo Fisher Scientific |
| 3 | `"Profiler Plus"` | Thermo Fisher Scientific |
| 4 | `"PowerPlex 18"` | Promega Corporation |
| 5 | `"PowerPlex 21"` | Promega Corporation |
| 6 | `"PowerPlex 16"` | Promega Corporation |
| 7 | `"DNATyper 15"` | Ministry of Public Security (MPS PRC) |
| 8 | `"Sinofiler"` | Thermo Fisher Scientific |
| 9 | `"SGM Plus"` | Thermo Fisher Scientific |
| 10 | `" ACGU 17+1"` *(leading space)* | Health BioMed (PRC) |
| 11 | `" GoldenEye 16A "` *(leading & trailing spaces)* | Peoplespot (PRC) |
| 12 | `" DNATyper 15 Plus "` *(leading & trailing spaces)* | MPS PRC |
| 13 | `" GoldenEye 16BT"` *(leading space)* | Peoplespot (PRC) |
| 14 | `" GoldenEye 16C"` *(leading space)* | Peoplespot (PRC) |
| 15 | `" GoldenEye 20A "` *(leading & trailing spaces)* | Peoplespot (PRC) |

---

## 5. Document Integrity and Validation Rules

A valid DDEM v0.9 document must satisfy all of the following criteria:

1. **Namespace Conformance**: Root element must declare `xmlns="urn:NDNADImportFile-schema"`.
2. **Global Uniqueness of Specimens**: Each `<SPECIMENID>` within a file must be unique (`UNIQUE_SPEC`).
3. **Locus Uniqueness within Specimen**: Within any single `<SPECIMEN>`, no `<LOCUSNAME>` may appear more than once (`UNIQUE_LOCI`).
4. **Field Length Restrictions**:
   - `SPECIMENID`: 1–24 characters.
   - `DESTINATIONLAB` / `SOURCELAB`: 1–10 characters.
   - `SUBMITBYUSERID` / `READINGBY`: 1–20 characters.
   - `ALLELEVALUE`: 1–10 characters.
   - `CASEID` / `BATCHID`: 0–32 characters.
   - `SPECIMENCOMMENT`: 0–255 characters.
5. **Mandatory Sequence Order**: All child elements must appear in the exact chronological sequence prescribed by the XSD `<xsd:sequence>` definitions.

---

## 6. Complete XML Example

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
        <SPECIMENCOMMENT>Reference blood standard</SPECIMENCOMMENT>
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
