---
layout: post
title:  "Breast cancer gene signature in R"
date:   2016-04-10 18:32:45 -0700
author: Ling Cheng
categories: 'Data Science'

---

```R
source("http://bioconductor.org/biocLite.R")
```

```R
biocLite("breastCancerNKI")
biocLite("genefu")
library(genefu)
library(breastCancerNKI)
library(Biobase)

library(ggplot2)
```

```R
data(sig.oncotypedx)
data(nkis)
data(nki)
nki.data <- exprs(nki)
nki.clinical=data.frame(pData(nki))
```

The dataset nki.data contains gene expression data, with expression profiling for ~25,000 genes from 337 patients. This is typical of gene expression data, with a lot more features than samples. Talk about curse of dimensionality. The second dataset is nki.clinical, with 21 features from 337 patients. Features include tumor size, age, tumor grade, er (estrogen receptor) positive or negative, spread to lymph node or not. This is where I benefited from being a biologist. It is hard to decipher what these abbreviated feature names mean without some knowledge of cancer biology.


```R
dim(nki.data)
```




<ol class=list-inline>
	<li>24481</li>
	<li>337</li>
</ol>





```R
dim(nki.clinical)
```




<ol class=list-inline>
	<li>337</li>
	<li>21</li>
</ol>





```R
head(nki.clinical)
```




<table>
<thead><tr><th></th><th scope=col>samplename</th><th scope=col>dataset</th><th scope=col>series</th><th scope=col>id</th><th scope=col>filename</th><th scope=col>size</th><th scope=col>age</th><th scope=col>er</th><th scope=col>grade</th><th scope=col>pgr</th><th scope=col>ellip.h</th><th scope=col>brca.mutation</th><th scope=col>e.dmfs</th><th scope=col>t.dmfs</th><th scope=col>node</th><th scope=col>t.rfs</th><th scope=col>e.rfs</th><th scope=col>treatment</th><th scope=col>tissue</th><th scope=col>t.os</th><th scope=col>e.os</th></tr></thead>
<tbody>
	<tr><th scope=row>NKI_4</th><td>NKI_4</td><td>NKI</td><td>NKI</td><td>4</td><td>NA</td><td>2</td><td>41</td><td>1</td><td>3</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>4747</td><td>0</td><td>4747</td><td>0</td><td>0</td><td>1</td><td>4744</td><td>0</td></tr>
	<tr><th scope=row>NKI_6</th><td>NKI_6</td><td>NKI</td><td>NKI</td><td>6</td><td>NA</td><td>1.3</td><td>49</td><td>1</td><td>2</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>4075</td><td>0</td><td>4075</td><td>0</td><td>0</td><td>1</td><td>4072</td><td>0</td></tr>
	<tr><th scope=row>NKI_7</th><td>NKI_7</td><td>NKI</td><td>NKI</td><td>7</td><td>NA</td><td>2</td><td>46</td><td>0</td><td>1</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>3703</td><td>0</td><td>3703</td><td>0</td><td>0</td><td>1</td><td>3700</td><td>0</td></tr>
	<tr><th scope=row>NKI_8</th><td>NKI_8</td><td>NKI</td><td>NKI</td><td>8</td><td>NA</td><td>2.8</td><td>48</td><td>0</td><td>3</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>3215</td><td>0</td><td>3215</td><td>0</td><td>0</td><td>1</td><td>3213</td><td>0</td></tr>
	<tr><th scope=row>NKI_9</th><td>NKI_9</td><td>NKI</td><td>NKI</td><td>9</td><td>NA</td><td>1.5</td><td>48</td><td>1</td><td>3</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>3760</td><td>0</td><td>3760</td><td>0</td><td>0</td><td>1</td><td>3757</td><td>0</td></tr>
	<tr><th scope=row>NKI_11</th><td>NKI_11</td><td>NKI</td><td>NKI</td><td>11</td><td>NA</td><td>2.2</td><td>37</td><td>1</td><td>3</td><td>NA</td><td>⋯</td><td>0</td><td>0</td><td>2120</td><td>0</td><td>2120</td><td>0</td><td>0</td><td>1</td><td>2119</td><td>0</td></tr>
</tbody>
</table>




Look through the features of nki.clinical. There are several possible target variables. 'dmfs' stands for distant metastasis free survival. 'rfs' stands for remission-free survial, and 'os' for overall survival. 't' for number of days, as continuous variable, and 'e' for events, as binary variable. These patients have been followed up for at least five years. After reading the original paper and other literature, e.dmfs is selected as target variable.


```R
#plot here but don't need it for the blog
nki.annot <- nki@featureData@data
rs.nki <- oncotypedx(data = t(nki.data), annot = subset(nki.annot, !is.na(EntrezGene.ID)), do.mapping = TRUE)
```


```R
nki.alldata <- subset(merge(nki.clinical, data.frame(samplename = names(rs.nki$score), oncotypedx = rs.nki$score), by = "samplename"), !is.na(oncotypedx))
```


```R
plot(survfit(Surv(t.rfs, e.rfs) ~ cut(oncotypedx, c(-1, 25, 50, 75, 101)), data = nki.alldata))
```


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_11_0.svg)


About 80% of all breast cancer are estrogen-receptor-positive. Tumors that are ER-positive are more likely to respond to hormone therapy and generally have better surival rates than ER-engative tumors. This dataset showed the same trend. ER-positive patients have significantly better outcomes, as shown by ttest and density plot


```R
nki.clinical$er = as.factor(nki.clinical$er)
t.test(t.dmfs~er, data=nki.clinical)
```





    	Welch Two Sample t-test

    data:  t.dmfs by er
    t = -3.5467, df = 118.15, p-value = 0.0005602
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -1113.2280  -315.5176
    sample estimates:
    mean in group 0 mean in group 1
           1997.203        2711.576





```R
ggplot(nki.clinical,aes(x=t.dmfs,fill=er))+geom_density(alpha=0.5)
```

    Warning message:
    : Removed 18 rows containing non-finite values (stat_density).


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_14_1.svg)


BRCA1 and BRCA2 are the best-known genes linked to breast cancer risk, and receive a lot of press coverage lately. This dataset has inforamtion on BRCA mutations for 117 patients, with 18 patients with BRCA1 mutation and 2 patients with BRCA2 mutation. Breast cancer linked with BRCA mutations are generally estrogen-receptor-negative and have less treatment options. Plotting of five-year dmfs showed strongly-reduced surivial.


```R
nki.clinical$brca.mutation=as.factor(nki.clinical$brca.mutation)
summary(nki.clinical$brca.mutation)
```




<dl class=dl-horizontal>
	<dt>0</dt>
		<dd>97</dd>
	<dt>1</dt>
		<dd>18</dd>
	<dt>2</dt>
		<dd>2</dd>
	<dt>NA's</dt>
		<dd>220</dd>
</dl>




Among 337 patients, 117 were screened for BRCA mutations. 97 had no BRCA mutations, 18 had BRCA1 mutation and 2 had BRCA2 mutation


```R
ggplot(nki.clinical,aes(x=t.dmfs))+geom_histogram()+facet_wrap(~brca.mutation)
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
    Warning message:
    : Removed 18 rows containing non-finite values (stat_bin).


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_18_1.svg)


Lymph node status is considered an important indicator of cancer status and survial. Survival of patients with negative lymph nodes is genereally better than those with positive lymph nodes. But in this study, t.dfms showed no significant difference between node-positive and node-negative patients. Gene expression signature, as discussed below, is a better indicator of cancer status and survival.


```R
nki.clinical$node=as.factor(nki.clinical$node)
t.test(t.dmfs~node,data=nki.clinical)
```





    	Welch Two Sample t-test

    data:  t.dmfs by node
    t = -0.42231, df = 311.31, p-value = 0.6731
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -408.6473  264.2286
    sample estimates:
    mean in group 0 mean in group 1
           2513.263        2585.472




Density plot shows a slight left shift for node-positive patients, consistent with the insignicant p-value.


```R
ggplot(nki.clinical,aes(x=t.dmfs,fill=node))+geom_density(alpha=0.5)
```

    Warning message:
    : Removed 18 rows containing non-finite values (stat_density).


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_22_1.svg)


Tumor grade is an important indicator of cancer survial. This dataset confirmed that based on anova test and density curve


```R
nki.clinical$grade = as.factor(nki.clinical$grade)
aov(t.dmfs~grade,data=nki.clinical)
```




    Call:
       aov(formula = t.dmfs ~ grade, data = nki.clinical)

    Terms:
                        grade Residuals
    Sum of Squares   26584290 714684554
    Deg. of Freedom         2       316

    Residual standard error: 1503.882
    Estimated effects may be unbalanced
    18 observations deleted due to missingness




```R
model.tables(aov(t.dmfs~grade,data=nki.clinical))
```




    Tables of effects

     grade
            1   2      3
        344.6 159 -331.2
    rep  78.0 108  133.0




```R
model.lm = lm(t.dmfs~grade,data=nki.clinical)
```


```R
anova(model.lm)
```




<table>
<thead><tr><th></th><th scope=col>Df</th><th scope=col>Sum Sq</th><th scope=col>Mean Sq</th><th scope=col>F value</th><th scope=col>Pr(>F)</th></tr></thead>
<tbody>
	<tr><th scope=row>grade</th><td>2</td><td>26584290</td><td>13292145</td><td>5.877163</td><td>0.003118219</td></tr>
	<tr><th scope=row>Residuals</th><td>316</td><td>714684554</td><td>2261660</td><td>NA</td><td>NA</td></tr>
</tbody>
</table>





```R
ggplot(nki.clinical,aes(x=t.dmfs,fill=grade))+geom_density(alpha=0.5)
```

    Warning message:
    : Removed 18 rows containing non-finite values (stat_density).


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_28_1.svg)



```R

```

Through supervised classification and leave-one-out sequential forward selection, van 't Veer et al. identified the 70-gene signature and acheivied >90% sensitivity with this signature. I used sig.gene70 for feature selection.


```R
data(sig.gene70)
dim(sig.gene70)
nki.gene70 = nki.data[sig.gene70$probe,]
dim(nki.gene70)
```




<ol class=list-inline>
	<li>70</li>
	<li>9</li>
</ol>







<ol class=list-inline>
	<li>70</li>
	<li>337</li>
</ol>





```R
head(nki.gene70)
```




<table>
<thead><tr><th></th><th scope=col>NKI_4</th><th scope=col>NKI_6</th><th scope=col>NKI_7</th><th scope=col>NKI_8</th><th scope=col>NKI_9</th><th scope=col>NKI_11</th><th scope=col>NKI_12</th><th scope=col>NKI_13</th><th scope=col>NKI_14</th><th scope=col>NKI_17</th><th scope=col></th><th scope=col>NKI_393</th><th scope=col>NKI_394</th><th scope=col>NKI_395</th><th scope=col>NKI_396</th><th scope=col>NKI_397</th><th scope=col>NKI_398</th><th scope=col>NKI_401</th><th scope=col>NKI_402</th><th scope=col>NKI_403</th><th scope=col>NKI_404</th></tr></thead>
<tbody>
	<tr><th scope=row>NM_003748</th><td>0.231 </td><td>0.159 </td><td>0.047 </td><td>-0.123</td><td>0.072 </td><td>0.011 </td><td>0.074 </td><td>0.324 </td><td>0.126 </td><td>0.126 </td><td>⋯     </td><td>0.138 </td><td>-0.105</td><td>0.026 </td><td>0.369 </td><td>-0.2  </td><td>-0.156</td><td>0.158 </td><td>-0.009</td><td>-0.175</td><td>-0.326</td></tr>
	<tr><th scope=row>NM_003862</th><td>-0.033</td><td>0.219 </td><td>0.77  </td><td>-0.418</td><td>0.186 </td><td>0.336 </td><td>-0.439</td><td>-0.095</td><td>0.226 </td><td>0.568 </td><td>⋯     </td><td>-0.066</td><td>-0.134</td><td>0.082 </td><td>-0.565</td><td>0.32  </td><td>-0.561</td><td>-0.521</td><td>-0.306</td><td>0.107 </td><td>-0.001</td></tr>
	<tr><th scope=row>Contig32125_RC</th><td>0.109 </td><td>0.352 </td><td>0.183 </td><td>-0.082</td><td>0.047 </td><td>-0.024</td><td>NaN   </td><td>0.376 </td><td>-0.032</td><td>0.312 </td><td>⋯     </td><td>-0.141</td><td>0.193 </td><td>0.018 </td><td>-0.142</td><td>-0.059</td><td>-0.252</td><td>0.109 </td><td>-0.064</td><td>-0.006</td><td>0.016 </td></tr>
	<tr><th scope=row>U82987</th><td>0.125 </td><td>0.394 </td><td>0.007 </td><td>0.139 </td><td>0.278 </td><td>-0.026</td><td>-0.235</td><td>0.179 </td><td>NaN   </td><td>0.235 </td><td>⋯     </td><td>0.081 </td><td>0.164 </td><td>0.164 </td><td>0.112 </td><td>0.111 </td><td>-0.184</td><td>-0.026</td><td>-0.105</td><td>0.114 </td><td>0.175 </td></tr>
	<tr><th scope=row>AB037863</th><td>0.456 </td><td>0.319 </td><td>0.323 </td><td>-0.253</td><td>0.294 </td><td>-0.303</td><td>-0.211</td><td>0.121 </td><td>-0.069</td><td>0.433 </td><td>⋯     </td><td>0.017 </td><td>0.319 </td><td>-0.048</td><td>-0.39 </td><td>-0.126</td><td>-0.171</td><td>0.179 </td><td>-0.297</td><td>-0.081</td><td>0.116 </td></tr>
	<tr><th scope=row>NM_020974</th><td>-0.628</td><td>0.703 </td><td>-0.218</td><td>-0.098</td><td>0.711 </td><td>-0.424</td><td>-1.18 </td><td>0.075 </td><td>0.881 </td><td>0.301 </td><td>⋯     </td><td>-0.66 </td><td>0.07  </td><td>0.091 </td><td>-1.201</td><td>0.044 </td><td>-0.934</td><td>0.963 </td><td>-0.846</td><td>0.464 </td><td>0.892 </td></tr>
</tbody>
</table>





```R
#transpose the matrix and remove NA
nki.gene70 = t(nki.gene70)
```


```R
nki.X = nki.gene70[complete.cases(nki.gene70),]
dim(nki.X)
```




<ol class=list-inline>
	<li>318</li>
	<li>70</li>
</ol>





```R
nki.Y = nki.clinical[complete.cases(nki.gene70),]
dim(nki.Y)
```




<ol class=list-inline>
	<li>318</li>
	<li>21</li>
</ol>





```R
#remove NA cases for e.dmfs
nki.Y2 = nki.Y[complete.cases(nki.Y$e.dmfs),]
dim(nki.Y2)
nki.X2 = nki.X[complete.cases(nki.Y$e.dmfs),]
dim(nki.X2)
```




<ol class=list-inline>
	<li>300</li>
	<li>21</li>
</ol>







<ol class=list-inline>
	<li>300</li>
	<li>70</li>
</ol>





```R
y=as.factor(nki.Y2$e.dmfs)
X=data.frame(nki.X2,y)
dim(X)
```




<ol class=list-inline>
	<li>300</li>
	<li>71</li>
</ol>







    NULL




```R
head(X)
```




<table>
<thead><tr><th></th><th scope=col>NM_003748</th><th scope=col>NM_003862</th><th scope=col>Contig32125_RC</th><th scope=col>U82987</th><th scope=col>AB037863</th><th scope=col>NM_020974</th><th scope=col>Contig55377_RC</th><th scope=col>NM_003882</th><th scope=col>NM_000849</th><th scope=col>Contig48328_RC</th><th scope=col>ellip.h</th><th scope=col>NM_020188</th><th scope=col>AL137718</th><th scope=col>Contig28552_RC</th><th scope=col>Contig38288_RC</th><th scope=col>AA555029_RC</th><th scope=col>NM_016359</th><th scope=col>Contig46218_RC</th><th scope=col>Contig63649_RC</th><th scope=col>AL080059</th><th scope=col>y</th></tr></thead>
<tbody>
	<tr><th scope=row>NKI_4</th><td>0.231</td><td>-0.033</td><td>0.109</td><td>0.125</td><td>0.456</td><td>-0.628</td><td>0.024</td><td>-0.276</td><td>-0.233</td><td>0.426</td><td>⋯</td><td>0.314</td><td>-0.038</td><td>-0.036</td><td>-0.037</td><td>-0.143</td><td>-0.01</td><td>-0.059</td><td>-0.077</td><td>0.235</td><td>0</td></tr>
	<tr><th scope=row>NKI_6</th><td>0.159</td><td>0.219</td><td>0.352</td><td>0.394</td><td>0.319</td><td>0.703</td><td>0.174</td><td>0.134</td><td>0.468</td><td>0.775</td><td>⋯</td><td>-0.279</td><td>-0.194</td><td>-0.365</td><td>0.078</td><td>-0.269</td><td>-0.349</td><td>-0.326</td><td>-0.147</td><td>-0.702</td><td>0</td></tr>
	<tr><th scope=row>NKI_7</th><td>0.047</td><td>0.77</td><td>0.183</td><td>0.007</td><td>0.323</td><td>-0.218</td><td>0.071</td><td>0.377</td><td>0.289</td><td>-0.243</td><td>⋯</td><td>-0.349</td><td>-0.134</td><td>-0.319</td><td>-0.092</td><td>-0.148</td><td>-0.452</td><td>-0.295</td><td>-0.27</td><td>-0.332</td><td>0</td></tr>
	<tr><th scope=row>NKI_8</th><td>-0.123</td><td>-0.418</td><td>-0.082</td><td>0.139</td><td>-0.253</td><td>-0.098</td><td>-0.143</td><td>-0.163</td><td>-0.383</td><td>-0.32</td><td>⋯</td><td>0.094</td><td>0.008</td><td>0.12</td><td>0.131</td><td>-0.14</td><td>0.006</td><td>-0.017</td><td>0.138</td><td>-0.234</td><td>0</td></tr>
	<tr><th scope=row>NKI_9</th><td>0.072</td><td>0.186</td><td>0.047</td><td>0.278</td><td>0.294</td><td>0.711</td><td>-0.168</td><td>0.195</td><td>0.501</td><td>0.322</td><td>⋯</td><td>-0.127</td><td>0.012</td><td>-0.042</td><td>-0.133</td><td>-0.18</td><td>-0.209</td><td>-0.118</td><td>0.12</td><td>-0.553</td><td>0</td></tr>
	<tr><th scope=row>NKI_11</th><td>0.011</td><td>0.336</td><td>-0.024</td><td>-0.026</td><td>-0.303</td><td>-0.424</td><td>-0.084</td><td>0.055</td><td>0.202</td><td>-0.541</td><td>⋯</td><td>0.125</td><td>-0.006</td><td>-0.047</td><td>-0.145</td><td>0.034</td><td>0.021</td><td>-0.065</td><td>-0.358</td><td>-0.351</td><td>0</td></tr>
</tbody>
</table>




First test randomforest model


```R
#install.packages('randomForest',"/home/user/anaconda3/lib/R/library")
library(randomForest)
```

    randomForest 4.6-12
    Type rfNews() to see new features/changes/bug fixes.

    Attaching package: ‘randomForest’

    The following object is masked from ‘package:ggplot2’:

        margin

    The following object is masked from ‘package:Biobase’:

        combine

    The following object is masked from ‘package:BiocGenerics’:

        combine




```R
train = sample (1: nrow(X), nrow(X)/2)
length(train)
```




150




```R
rf.nki = randomForest(y~.,data=X,subset=train,mtyr=9)
```


```R
rf.pred = predict(rf.nki, newdata=X[-train,])
```


```R
table(rf.pred,y[-train])
```





    rf.pred  0  1
          0 85 42
          1  9 14




```R
mean(rf.pred==y[-train])
```




0.66




```R
#sesitivity/recall
14/(14+42)
```




0.25



SVM model


```R
library(e1071)
```


```R
#use the entire dataset and cross-validation
tune1=tune(svm,y~.,data=X,kernel='linear',ranges=list(cost=c(0.0001,0.001,0.01,0.1,1,10)))
```


```R
summary(tune1)
```





    Parameter tuning of ‘svm’:

    - sampling method: 10-fold cross validation

    - best parameters:
     cost
     0.01

    - best performance: 0.3

    - Detailed performance results:
       cost     error dispersion
    1 1e-04 0.3400000 0.04097575
    2 1e-03 0.3400000 0.04097575
    3 1e-02 0.3000000 0.08606630
    4 1e-01 0.3100000 0.08172522
    5 1e+00 0.3166667 0.08202679
    6 1e+01 0.3300000 0.10116604





```R
#use c=0.01 for train/test split
svm1=svm(y~.,data=X[train,],kernel='linear',cost=0.01)
```


```R
svm1.pred = predict(svm1, newdata=X[-train,])
```


```R
table(svm1.pred,y[-train])
```





    svm1.pred  0  1
            0 90 49
            1  4  7




```R
mean(svm1.pred==y[-train])
```




0.646666666666667




```R
#sensitivity/reacall
7/(7+49)
```




0.125




```R
svm2=svm(y~.,data=X[train,],kernel='linear',cost=0.1)
svm2.pred = predict(svm2, newdata=X[-train,])
table(svm2.pred,y[-train])
```





    svm2.pred  0  1
            0 70 30
            1 24 26




```R
mean(svm2.pred==y[-train])
```




0.64




```R
#sensitivity/recall
26/(26+30)
```




0.464285714285714




```R
#polynomial kernels
tune2=tune(svm,y~.,data=X,kernel='polynomial',ranges=list(degree=c(1,2,3),cost=c(0.01,0.1,1)))
```


```R
summary(tune2)
```





    Parameter tuning of ‘svm’:

    - sampling method: 10-fold cross validation

    - best parameters:
     degree cost
          1    1

    - best performance: 0.31

    - Detailed performance results:
      degree cost     error dispersion
    1      1 0.01 0.3400000 0.07665861
    2      2 0.01 0.3400000 0.07665861
    3      3 0.01 0.3400000 0.07665861
    4      1 0.10 0.3400000 0.07665861
    5      2 0.10 0.3400000 0.07665861
    6      3 0.10 0.3400000 0.07665861
    7      1 1.00 0.3100000 0.06488356
    8      2 1.00 0.3433333 0.08020037
    9      3 1.00 0.3333333 0.08748898





```R
tune3 = tune(svm,y~.,data=X,kernel='radial',ranges=list(cost=c(0.01,0.1,1,10,100),gamma=c(0.5,1,2,3,4)))
```


```R
summary(tune3)
```





    Parameter tuning of ‘svm’:

    - sampling method: 10-fold cross validation

    - best parameters:
     cost gamma
     0.01   0.5

    - best performance: 0.34

    - Detailed performance results:
        cost gamma error dispersion
    1  1e-02   0.5  0.34 0.07503086
    2  1e-01   0.5  0.34 0.07503086
    3  1e+00   0.5  0.34 0.07503086
    4  1e+01   0.5  0.34 0.07503086
    5  1e+02   0.5  0.34 0.07503086
    6  1e-02   1.0  0.34 0.07503086
    7  1e-01   1.0  0.34 0.07503086
    8  1e+00   1.0  0.34 0.07503086
    9  1e+01   1.0  0.34 0.07503086
    10 1e+02   1.0  0.34 0.07503086
    11 1e-02   2.0  0.34 0.07503086
    12 1e-01   2.0  0.34 0.07503086
    13 1e+00   2.0  0.34 0.07503086
    14 1e+01   2.0  0.34 0.07503086
    15 1e+02   2.0  0.34 0.07503086
    16 1e-02   3.0  0.34 0.07503086
    17 1e-01   3.0  0.34 0.07503086
    18 1e+00   3.0  0.34 0.07503086
    19 1e+01   3.0  0.34 0.07503086
    20 1e+02   3.0  0.34 0.07503086
    21 1e-02   4.0  0.34 0.07503086
    22 1e-01   4.0  0.34 0.07503086
    23 1e+00   4.0  0.34 0.07503086
    24 1e+01   4.0  0.34 0.07503086
    25 1e+02   4.0  0.34 0.07503086





```R
library(ROCR)
svm.nki =svm(y~.,data=X,kernel='linear',cost=0.1,probability=TRUE)
pred.svm =predict(svm.nki,X,probability=TRUE)
pred.svm.roc=prediction(attr(pred.svm,'probabilities')[,2],y)
perf.svm=performance(pred.svm.roc,'tpr','fpr')

rf.nki =randomForest(y~.,data=X,mtyr=9)
pred.rf.roc=prediction(as.vector(rf.nki$votes[,2]),y)
perf.rf=performance(pred.rf.roc,'tpr','fpr')
```

    Loading required package: gplots
    Warning message:
    : package ‘gplots’ was built under R version 3.2.4
    Attaching package: ‘gplots’

    The following object is masked from ‘package:stats’:

        lowess




```R
plot(perf.svm,col=1)
plot(perf.rf,add=TRUE,col=2)
legend('right',c('svm','randomForest'),text.col=c(1,2))
```


![svg](NKIbreastcancer_1_files/NKIbreastcancer_1_64_0.svg)



```R

```


```R

```


```R

```
