# MLCFforTracking

## Multiple Local Correlation Filters for Robust Visual Tracking

Correlation filter based trackers have been successfully applied to visual object tracking. These methods utilize a periodic assumption of the training samples that also introduces unwanted boundary effects. Recently, Spatially Regularized Correlation Filters (SRDCF) solved this issue by introducing a penalization term to the filter coefficients, enabling the boundary effects to be effectively reduced by assigning higher weights to the background. However, since the scale of the tracking target is variable, it is difficult to design the penalize ratio for the filter coefficients, resulting in poor robustness when the object is partially occluded. In this paper, we investigated SRDCF and present a novel region-based tracking algorithm using Multiple Local Correlation Filters (MLCF). In our tracking framework, all correlation filters are learned over the training samples under a circulant structure, and each correlation filter is associated with a different desired confidence map to achieve region-based tracking. 

**Our main contributions include the following:**
- A novel region-based tracking approach using Multiple Local Correlation Filters (MLCF). The MLCF model tracks each separate local region by learning multiple local correlation filters without cropping the object into parts.
- Incorporating region-based tracking strategy into the local correlation filters. Our tracker learns more reliable features from the target area. 
- Using PSR to adaptively update the target to avoid filters being contaminated when the tracker loses the object.

Comprehensive experiments were conducted on two large-scale benchmark datasets: OTB-2015 and OTB-2013, and the experimental results demonstrate that the proposed MLCF tracker performs competitively when compared with state-of-the-art methods.

## Comparisons
**qualitative**
Comparisons of our approach with state-of-the-art trackers in challenging situations of rotations (#1026 to #1078) and partial occlusions (#321 to #372) on the Lemming sequence [1]. Our tracker takes advantage of a region-based tracking strategy and correlation filters; it performs superiorly when compared with other trackers.
<img src="https://github.com/djzgroup/MLCFforTracking/blob/master/img/Comparisons.jpg" width="800">

**quantitative**
We do comprehensive comparison between our approach and the baseline tracker SRDCF. The following figure shows the DP and OP results of MLCF and SRDCF on OTB-2015 respectively, for clarity, we reported 11 attribute-based evaluation results. (Comparison with baseline tracker SRDCF on OTB-2015 benchmark sequences using OP at a threshold of 0.5, DP at a threshold of 20 pixels.)
<img src="https://github.com/djzgroup/MLCFforTracking/blob/master/img/baseline.jpg" width="700">

## CODE
The code is developed in the matlab environment.
```bash
You can verify the method of this article by executing “run_OUR.m”.
```
## References
- [1] Wu Y, Lim J, Yang M H. Object Tracking Benchmark[J]. IEEE transactions on pattern analysis and machine intelligence, 2015, 37(9): 1834.

## Acknowledgment
This work was supported by the National Natural Science Foundation of China under Grant 61702350 and Grant 61802355.
