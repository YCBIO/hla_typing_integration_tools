### Description
This tool is used to integrate the typing results of OptiType and Polysolver

### Function
Although OptiType has high typing accuracy for control samples, however, when we used the OptiType for HLA typing, some patients had inconsistent typing results in tumor and control samples. We utilize the Polysolver for a secondary analysis of these inconsistent samples. All inconsistent typing can find consistent typing in Polysolver analysis results, part of the tumor sample typing results of OptiType are consistent with Polysolver typing results.
Therefore, we developed a scoring algorithm to integrate the OptiType and Polysolver typing results of tumor and control samples to obtain the most reliable HLA typing. We will verify the integration results of our scoring program in the future. 

### Usage：

```bash
python3 Hla_typing_integration.py \
--on  /mnt/cfs/project/test_freshman/yijian/pipline/hla/test_data/P36-N_N_WES_P36-N_result.tsv  \
--ot  /mnt/cfs/project/test_freshman/yijian/pipline/hla/test_data/P36-T_T_WES_P36-T_result.tsv \
--pn /mnt/cfs/project/test_freshman/yijian/pipline/hla/test_data/P36-N.winners.hla.nofreq.txt \
--pt  /mnt/cfs/project/test_freshman/yijian/pipline/hla/test_data/P36-T.winners.hla.nofreq.txt \
-r /mnt/cfs/project/test_freshman/yijian/pipline/hla/test_data/P36.intgration.hla.txt
```
 
### parameter：<br />
  - -on, required, optitype typing result for normal sample
  - -ot, required, optitype typing result for tumor sample
  - -pn, required, polysolver typing result for normal sample 
  - -pt, required, polysolver typing result for tumor sample
  - -r, required, integrate result file，line 1 is the result after integration, and line 2-5 is the original typing result.

