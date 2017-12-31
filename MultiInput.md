
**We have seen that graphical models provide a compact way to define joint probability distribu- tions. Given such a joint distribution, what can we do with it? The main use for such a joint distribution is to perform probabilistic inference. This refers to the task of estimating unknown quantities from known quantities. **


### sigmoid and non-linearity

when the distribution differs, QQ plot is not linear


requirements:
1. transformed value range from 0 to 1
2. trasform is monotonic, keep rankings unchanged
3. non-linear

one way is to use a threshold function, with multiple thresholds to perserve as much information as possible, taking a extreme situation, when every single point is associated with a threshold, it's the ranking(or cdf) that point.

when using ${cdf}_x$ as transform function, the transformed values ${cdf}_x(x)$ are uniform distributed, satisfy requirements 1 and 2, but can't handle non-linearity.


Sigmoid functions are often used in artificial neural networks to
introduce nonlinearity in the model. 

The purpose of the activation function is to introduce non-linearity into the network
in turn, this allows you to model a response variable (aka target variable, class label, or score) that varies non-linearly with its explanatory variables

non-linear means that the output cannot be reproduced from a linear combination of the inputs (which is not the same as output that renders to a straight line--the word for this is affine).

$y = \frac{x-\mu_x}{\sigma_x}$ shows a linear relation between x and y

WheninfluencecanflowfromXtoY viaZ,wesaythatthetrailX⌦Z⌦Y isactive. The results of our analysis for active two-edge trails are summarized thus:

CausaltrailX!Z!Y:activeifandonlyifZisnotobserved.
• Evidential trail X Z Y : active if and only if Z is not observed.
• Common cause X Z ! Y : active if and only if Z is not observed.
• Common e ect X ! Z Y : active if and only if either Z or one of Z’s descendants is observed.


## Chapter 5 is useful in constructing graphical models, but it's not current focus