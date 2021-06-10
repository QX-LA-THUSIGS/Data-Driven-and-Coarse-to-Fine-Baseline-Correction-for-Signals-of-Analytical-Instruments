# Data-Driven-and-Coarse-to-Fine-Baseline-Correction-for-Signals-of-Analytical-Instruments
## You can use this program to realize the baseline adaptive correction of signals from a variety of analytical instruments,  including but not limited to mass spectrometers, ion mobility spectrometers, and chromatographs.  

__The algorithm overcomes the mode-mixing problem of empirical mode decomposition algorithm by adaptively locating and removing high-amplitude spectral peaks. Through qualitative and quantitative analysis, compared with the traditional least squares fitting and sparse representation, the algorithm based on empirical mode decomposition (DD-CF) has better advantages in processing time and baseline fitting effect. At the same time, the biggest feature of this algorithm is that it can realize the data-driven baseline correction of mass spectrometer, chromatograph and ion mobility spectrum without user intervention. Compared with traditional algorithms, DD-CF algorithm has stronger adaptive ability.__  

## Attention
  
__This algorithm uses MATLAB's built-in function EMD, which is only supported in version 2018a and later versions. Therefore, we recommend using a newer version of MATLAB. If you can only use MATLAB before 2018a, we recommend using the method mentioned in [this article](https://www.programmersought.com/article/94925533430/) to download and install the independent EMD toolkit.There are many similar EMD toolkit installation methods, you can also find and install it yourself. However, it should be noted that the calling method of the custom EMD function and the DD-CF algorithm mentioned in this article may be different. Please make corresponding improvements on the original algorithm according to the toolkit prompts.__



__For details of the algorithm, please refer to the original reference of this algorithm.(DOI:10.1016/j.aca.2021.338386)__  


## Algorithm introduction

![DD-CF algorithm flow chart](https://github.com/QX-LA-THUSIGS/Data-Driven-and-Coarse-to-Fine-Baseline-Correction-for-Signals-of-Analytical-Instruments/blob/main/picture/flow_chart.png)  
            _Figure1:DD-CF flow chart_

__step1:__  

  Polynomial fitting method to obtain rough baseline  


__step2:__ 

  Cubic spline interpolation to remove spectral peaks  
  
__step3:__  

  EMD and separate the baseline-led IMF  

__step4:__   

  Signal reconstruction  


## How to use  
1. The code folder contains two MATLAB script files：  
    DD_CF_v1.m is the function script of the DD-CF algorithm；  
    testCode.m is a test program that calls the DD_CF_v1.m script to verify the baseline correction effect of the DD-CF algorithm on the data in the folder data.  
2. The data folder contains multiple sets of real experimental data of mass spectra, ion mobility spectra, and chromatographs to verify the baseline correction effect of the algorithm. More mass spectrometry experimental data can be downloaded from Baidu Netdisk: [More mass spectrum data link](https://pan.baidu.com/s/1MQ5bopY8lhOMwYQmR9g_-w). Extraction code：f09i. 
3. You can use the DD_CF_v1 function __"[DBSig,baseline] = DD_CF_v1(Sig)"__ to reslize the baseline adaptive correction. where, DBSig is the signal which represents signal after baseline correction;baseline represents the fitting baseline;Sig represents signal needed to be removed baseline.

## Real chemical signal baseline correction example  
__1. Mass spectrum__  
![mass spectrum signal](https://github.com/QX-LA-THUSIGS/Data-Driven-and-Coarse-to-Fine-Baseline-Correction-for-Signals-of-Analytical-Instruments/blob/main/picture/ms.png)    
            _Figure2:the baseline correction of mass spectrum signal_  

__2. Chromatographic instrument__
![Chromatographic instrument signal](https://github.com/QX-LA-THUSIGS/Data-Driven-and-Coarse-to-Fine-Baseline-Correction-for-Signals-of-Analytical-Instruments/blob/main/picture/chromatographic.png)    
            _Figure3:the baseline correction of chromatographic instrument signal_  

__3. Mobility spectrum__
![Mobility spectrum signal](https://github.com/QX-LA-THUSIGS/Data-Driven-and-Coarse-to-Fine-Baseline-Correction-for-Signals-of-Analytical-Instruments/blob/main/picture/migration_spectrum.png)    
            _Figure4:the baseline correction of mobility spectrum signal_ 


## Statement

The copyright belongs to the Analytical Instrument Group of Shenzhen International Graduate School of Tsinghua University.  
This program is designed and maintained by the QX-LA-THUSIGS team.  

__If you need to cite or improve this algorithm, please declare the original document of this algorithm in your work.)__  
__References:__  
Xiangchun Xu, et.al, Data-Driven and Coarse-to-Fine Baseline Correction for Signals of Analytical Instruments[J],Analytica Chimica Acta,2021,1157(338386):1873-4324.（DOI:10.1016/j.aca.2021.338386）
