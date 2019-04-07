# MLCFforTracking

## Multiple Local Correlation Filters for Robust Visual Tracking

Correlation filter based trackers have been successfully applied to visual object tracking. These methods utilize a periodic assumption of the training samples that also introduces unwanted boundary effects. Recently, Spatially Regularized Correlation Filters (SRDCF) solved this issue by introducing a penalization term to the filter coefficients, enabling the boundary effects to be effectively reduced by assigning higher weights to the background. However, since the scale of the tracking target is variable, it is difficult to design the penalize ratio for the filter coefficients, resulting in poor robustness when the object is partially occluded. In this paper, we investigated SRDCF and present a novel region-based tracking algorithm using Multiple Local Correlation Filters (MLCF). In our tracking framework, all correlation filters are learned over the training samples under a circulant structure, and each correlation filter is associated with a different desired confidence map to achieve region-based tracking. 

**Our main contributions include the following:**
- A novel region-based tracking approach using Multiple Local Correlation Filters (MLCF). The MLCF model tracks each separate local region by learning multiple local correlation filters without cropping the object into parts.
- Incorporating region-based tracking strategy into the local correlation filters. Our tracker learns more reliable features from the target area. 
- Using PSR to adaptively update the target to avoid filters being contaminated when the tracker loses the object.

Comprehensive experiments were conducted on two large-scale benchmark datasets: OTB-2015 and OTB-2013, and the experimental results demonstrate that the proposed MLCF tracker performs competitively when compared with state-of-the-art methods.

## CODE
The code is developed in the matlab environment.
```bash
You can verify the method of this article by executing “test-ucf24.py”.
```

## Acknowledgment
This work was supported by the National Natural Science Foundation of China under Grant 61702350 and Grant 61802355.
