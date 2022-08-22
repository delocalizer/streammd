# Goal

Compare GS results using streammd vs picard MD.

# Method

Updated tasks in `genomeinfo/production/wdl/branches/crl_streammd` to use
qname sorted chunks and streammd as a duplicate marker.

Ran the workflow through cromwell.

# Results

## Target

The GS_38_maf regression test targets are:
```
SNP_GP_expected=61654
SNP_GPC_expected=369
SNP_SP_expected=33707
SNP_SPC_expected=389
INDEL_GP_expected=4212
INDEL_GPC_expected=110
INDEL_SP_expected=144
```

The counts we get from the streammd workflow:
```
gp.snpMaf.maf=61596
gpc.snpMaf.maf=369
sp.snpMaf.maf=33690
spc.snpMaf.maf=389
gp.indelMaf.maf=4212
gpc.indelMaf.maf=110
sp.indelMaf.maf=144
```

So pretty close by line count; how about by exact position match?
```
grep -f <(diff <(grep -v '^#' Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.gp.snpMaf.maf|cut -f5,6,7) <(grep -v '^#' /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.gp.snpMaf.maf |cut -f5,6,7)|grep  '^>'|sed 's/^> //') Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.snpMaf.maf |wc -l
105
[conradL@hpcnode064 execution]$ grep -f <(diff <(grep -v '^#' Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.gp.snpMaf.maf|cut -f5,6,7) <(grep -v '^#' /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.gp.snpMaf.maf |cut -f5,6,7)|grep  '^<'|sed 's/^< //') /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.snpMaf.maf |wc -l
47

```
i.e. 105 positions (out of 61k approx) are in the picard wf but not streammd, 
and 47 positions are in the streammd but not picard.

What are the details?
```
grep -f <(diff <(grep -v '^#' Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.gp.snpMaf.maf|cut -f5,6,7) <(grep -v '^#' /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.gp.snpMaf.maf |cut -f5,6,7)|grep  '^>'|sed 's/^> //') Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.snpMaf.maf |cut -f 35|sort|uniq -c|sed 's/PASS,\|,PASS//g'
      9 5BP=1
      3 5BP=2
      1 5BP=3
      1 5BP=9
      5 COV,COV
      1 COV,COV,COV
      1 COV,COV
     11 COV
      1 COV;5BP=1
      1 COV;MR
      1 COV;SBIASCOV
      1 MIN,COV,COV
      3 MIN
      1 MIUN,COV,COV
      1 MR,COV,COV
      1 MR,MR
     15 MR
      1 MR;5BP=1
      5 NNS;MR
      1 5BP=1,MR
      1 5BP=1
      2 COV,COV,COV
      3 MR
      1 NNS;MR
     13 COV,COV
      1 COV,COV;MR
      3 MR
      1 MR
      1 SBIASALT,MR
      1 SBIASALT,COV,COV
      6 SBIASALT
      1 SBIASALT,SBIASALT;5BP=15
      1 SBIASALT;5BP=1
      1 SBIASALT;5BP=2,SBIASALT;5BP=5
      1 SBIASALT;MR,MR
      1 SBIASALT;NNS;MR;5BP=1
      1 SBIASCOV,SBIASCOV
      2 SBIASCOV;5BP=1

grep -f <(diff <(grep -v '^#' Colo-829.Colo-829-A.95a82ee1-50a4-4a64-8f00-f13e0ecf6039_vs_Colo-829.Colo-829-BL-B.dbb89b87-38f4-42f9-9687-9bf7dc1c2964.gp.snpMaf.maf|cut -f5,6,7) <(grep -v '^#' /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.gp.snpMaf.maf |cut -f5,6,7)|grep  '^<'|sed 's/^< //') /working/genomeinfo/cromwell-test/cromwell-executions/somaticDnaFastqToMaf/8d39dec5-6504-4a84-b027-84756ec24279/call-mergeMafs/execution/Colo-829.Colo-829-A.87b7f4af-dd45-4635-bb48-61cebeeb9429_vs_Colo-829.Colo-829-BL-B.e727f09f-8725-45c1-9009-adfde0be4b9f.snpMaf.maf |cut -f 35|sort|uniq -c|sed 's/,PASS\|PASS,//g'
      2 5BP=1
      1 5BP=1,SBIASALT;5BP=1
      1 5BP=2,COV,COV
      1 5BP=2
      1 5BP=2,SBIASALT;5BP=6
      1 COV,COV
      2 COV
      2 MIN
      9 MR
      1 MR;5BP=1
      1 MR;5BP=3
      2 NNS
      3 NNS;MR
      1 NNS;MR;5BP=1
      1 5BP=1
      1 MR
      6 COV,COV
      1 MR
      1 MR
      1 SBIASALT
      3 SBIASALT
      1 SBIASALT;5BP=2
      1 SBIASALT;MR;5BP=1
      2 SBIASCOV
      1 SBIASCOV;5BP=1
```
