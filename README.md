# Modelling Market Response for Axe Deodorants in Dutch Supermarkets

In this project, a real-life data set that contains data for several brands of deosprays for a number of supermarket chains in the Netherlands is used. Weekly sales data are available for eight brands (Dove, Fa, Nivea Rexona, Sanex, Vogue, 8x4, and Axe) across five well-known Dutch supermarket formulas (Albert Heijn, Edah, Super de Boer, Jumbo, C1000). The deospray brands differ in terms of their price positioning, which might be related to the tiers the brands are in. The supermarket formulas are chosen based on size/importance in the market, the positioning of the formula, but also based on data availability.

The detailed analysis and results can be found in the report. Codes are also provided. Below is the brief summary of the project.

### Description of the Dataset
For each of the eight brands and for each of the five supermarket formulas, 124 weekly data points are available (running from week 46 in 2003 until week 12 in 2006). For each brand / formula combination, the variables are:
  - Sales (volume);
  - Prices (euros per l00 ml) (as observed);
  - Regular price (euros per 100 ml) (proxy for the price as if there were no discounts);
  - Variable for ‘feature-only’: weighted distribution figure at the chain level for a ‘feature-only’ (no display, see below) in a certain week.
  - Variable for display-only: weighted distribution figure for display-only (no feature) in a certain week. 
  - Variable for combined feature and display: again a weighted distribution figure. It is measured as a dummy at the store level: 1 for combined use of feature and display, else 0.

Hence, at the store-level, we have four situations: (1) no promotion (feature-only, displayonly, feature & display all zero), (2) feature-only promotion: only this variable is 1, others
are zero, (3) display-only promotion: only this variable is 1, the others zero, and (4)
combined feature-and-display promotion: this variable is one, the others zero.

## PART A - Exploration, Specification, Estimation and Validation
In order to have a better understanding of the data, certain questions are answered regarding market share, retailers, competition, trends and seasonality, and other issues.

#### Data Analysis
Quarterly market share of Axe deodorants is extracted on general level and on retailer level. Sales of AXE deodorants are not distributed equally across the retailers. Then, price levels of AXE deodorants and of the competitors are analyzed. 2005 is important, because that is when the price war began. AXE deodorants dropped the price to a more stable level around 2.8€. Similar drop can be seen from the figures of other competitors. Price levels of AXE deodorants follow a similar pattern across the retailers. Therefore, price levels in different supermarkets follow the general trend. 

When promotions are analyzed, it is seen that all three kinds of promotion has asignificantly positive impact on sales for Albert Heijn and C-1000. These two supermarket chains constitute about 68% of AXE sales. Therefore, it is recommended to continue with both kind of promotions for these two supermarket chains to increase sales.

#### Response Model for Market Share
After the data analysis process, a response model is specified for the weekly market share of AXE deodorants in Albert Heijn supermarkets. Albert Heijn is selected, because it provides the most sales volume among the five supermarket chains in the data. It is assumed that the market is a duopoly – AXE and the others. Therefore, sales and price-related variables are averaged, and promotionrelated variables are summed.

With the model-ready data, three models are specified. One is linear, one is multiplicative and one is constrained. Constraint is necessary, because the decision variable (market share) has an upper bound of 1 (100%) and it cannot be more than that. After comparing these models, the linear model is found to be the best performer. The linear model can explain 78.39% of the variability in the weekly market share of Axe deodorants at Albert Heijn chains.

## PART B - Dynamic Response Models
This time, dynamics are included into the model, and the temporal nature of the competitive interactions is assessed, in terms of the promotional support that brands receive.

New wariables (relative price, dummies for quarters, price war etc.) are added and the model is modified. Now, the explainability increased to 79.92%. It is found that relative prices are a better indicator than absolute prices.

#### Lags
Lags are incorporated to the model to accommodate dynamic effects of the market. Only relative prices (2 lags) and combined promotions (1 lag) are found to have significant lagged effects. Then, with this model with lagged variables on top of the previous model, the explainability increased to 85.09%.

#### Most Effective Promotional Instruments
In general, the combined promotion levels for Axe deodorants have higher impact on the market share than relative price, meaning that the elasticity of the combined promotion instrument is higher.

#### Static vs. Dynamic
In terms of overall predictive performance, dynamic model yielded slightly better results, because the RMSE and MAE values of the dynamic model are smaller than those of the static model. In terms of explained variability, static model performed slightly better, because the R2 of the static model is greater than that of the dynamic model. However, the dynamic model captures the lagged effects for relative price and combined promotion levels for Axe deodorants. Therefore, the conclusion is that the dynamic model should be preferred over the static model, due to higher explainability for managerial decisions.

## PART C - Conclusions
The key insights are summarized below:
  - Even though the product is the same (deodorants), the differences in branding position for companies and supermarket chains, according to price and promotion levels, yield different outcomes in terms of sales volume.
  - For Axe deodorants, around 80% of the sales come from only 3 supermarket chains: Albert Heijn, C-1000 and Super De Boer. They provide higher service than some other chains. Therefore, future actions should mainly be applicable to high-service chains.
  - According to the prediction model, decreasing prices have positive impact on market share. Decreasing the prices down to the level of 2nd most expensive brand – or maybe even slightly
lower – might be beneficial in terms of sales and market share.
  - It can be seen from the prediction model coefficients that combined promotions (feature and display together) have much more impact on market share than price levels. Therefore, focusing more on promotions than on price might be beneficial, especially when and where competitors have less promotions.
  - The lagged effects for relative price and combined promotion levels for Axe deodorants
have significant impact on the market share. Therefore, dynamic model should be preferred
for forecasting purposes and for managerial decisions.
  - The elasticity of combined promotion levels is higher than that of relative price levels.
Therefore, it can be concluded that combined promotion level is a stronger instrument for
increasing the market share.
  - The impact of combined promotion levels, the elasticity, has a decaying effect, that is, the
effect of combined promotion at a given time affects less and less for the future market
share.
