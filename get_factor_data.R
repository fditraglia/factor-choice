#-----------------------------------------------------
# Liquidity Factor (liq_data_1962_2014.txt)
# Source: Lubos Pastor's Website
# (faculty.chicagobooth.edu/lubos.pastor/research/)
#-----------------------------------------------------
# DOCUMENTATION:
#
#   Liquidity Factor of Pastor & Stambaugh (JPE 2003)
#   Updated Through: December 2014
#  			
#   Column 1: Month	
#   Column 2: Levels of aggregate liquidity
#               (Figure 1; equation (5))		
#   Column 3: Innovations in aggregate liquidity 
#               (non-traded liquidity factor; 
#                equation (8); the main series)			
#   Column 4: Traded liquidity factor 
#               (LIQ_V, 10-1 portfolio return)			
# 			
#   NOTE: The traded factor is the value-weighted return
#         on the 10-1 portfolio from a sort on historical
#         liquidity betas. This procedure is simpler than
#         sorting on predicted betas (as in the original
#         study), and through 2014 it is similarly
#         successful at creating a spread in post-ranking 
#         betas. The traded	factor has a positive and
#         significant alpha through 2014, consistent with
#         liquidity risk being priced.			
#-----------------------------------------------------
liquidity_url <- "http://faculty.chicagobooth.edu/lubos.pastor/research/liq_data_1962_2014.txt"
liquidity <- read.delim(liquidity_url, comment.char = "%")
names(liquidity) <- c("month", "LIQ.agg", 
                      "LIQ.innov", "LIQ.v")
rm(liquidity_url)


#-----------------------------------------------------
# Momentum Factor (F-F_Momentum_Factor.zip)
# Source: Ken French's Website
# (mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp)
#-----------------------------------------------------
# DOCUMENTATION:
# 
#     Monthly Returns:  January 1927 - March 2015
#
#     We use six value-weight portfolios formed on size 
#     and prior (2-12) returns to construct Mom. The
#     portfolios, which are formed monthly, are the
#     intersections of 2 portfolios formed on size (market
#     equity, ME) and 3 portfolios formed on prior (2-12)
#     return. The monthly size breakpoint is the median
#     NYSE market equity. The monthly prior (2-12) return
#     breakpoints are the 30th and 70th NYSE percentiles.
#
#     Mom is the average return on the two high prior
#     return portfolios minus the average return on the 
#     two low prior return portfolios,
#
#                   Mom = 	1/2 (Small High + Big High)
#                         - 1/2(Small Low + Big Low). 	 
#
#     The six portfolios used to construct Mom each month
#     include NYSE, AMEX, and NASDAQ stocks with prior
#     return data. To be included in a portfolio for month
#     t (formed at the end of month t-1), a stock must have
#     a price for the end of month t-13 and a good return
#     for t-2. In addition, any missing returns from t-12
#     to t-3 must be -99.0, CRSP's code for a missing price.
#     Each included stock also must have ME for the end of
#     month t-1.
#-----------------------------------------------------
momentum_url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor.zip"
temp <- tempfile()
download.file(momentum_url, temp)
momentum <- read.table(unz(temp, "F-F_Momentum_Factor.TXT"),
                      na.strings = c("-999", "-99.99"), 
                      skip = 13, nrows = 1059) 
unlink(temp)
rm(momentum_url, temp)
momentum <- data.frame(month = rownames(momentum),
                       momentum)
rownames(momentum) <- NULL
names(momentum)[2] <- "MOM"


#-----------------------------------------------------
# Fama/French 5 Factors (2x3)
# Source: Ken French's Website
# (mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp)
#-----------------------------------------------------
# DOCUMENTATION:
#
#   Monthly Returns: July 1963 - March 2015
#
#   The Fama/French 5 factors (2x3) are constructed
#   using the 6 value-weight portfolios formed on size 
#   and book-to-market, the 6 value-weight portfolios 
#   formed on size and operating profitability, and 
#   the 6 value-weight portfolios formed on size and 
#   investment. (See the description of the 6 
#   size/book-to-market, size/operating profitability, 
#   size/investment portfolios.)
#
#   SMB (Small Minus Big) is the average return on the 
#   nine small stock portfolios minus the average 
#   return on the nine big stock portfolios,
#
#     SMB(B/M) = 1/3(Small Value + Small Neutral + Small Growth)
#              - 1/3(Big Value + Big Neutral + Big Growth).
#
#     SMB(OP) = 1/3(Small Robust + Small Neutral + Small Weak)
#             - 1/3(Big Robust + Big Neutral + Big Weak).
#
#     SMB(INV) = 1/3 (Small Conservative + Small Neutral + 
#                     Small Aggressive) - 1/3 (Big Conservative
#                     + Big Neutral + Big Aggressive).
#
#   SMB = 1/3 ( SMB(B/M) + SMB(OP) + SMB(INV) ).
#
#   HML (High Minus Low) is the average return on the 
#   two value portfolios minus the average return on 
#   the two growth portfolios,
#
#     HML = 1/2 (Small Value + Big Value)
#         - 1/2 (Small Growth + Big Growth). 	 
#
#   RMW (Robust Minus Weak) is the average return on 
#   the two robust operating profitability portfolios 
#   minus the average return on the two weak operating 
#   profitability portfolios,
#
#     RMW = 1/2 (Small Robust + Big Robust)
#         - 1/2 (Small Weak + Big Weak). 	 
#
#   CMA (Conservative Minus Aggressive) is the average 
#   return on the two conservative investment portfolios 
#   minus the average return on the two aggressive 
#   investment portfolios,
#
#     CMA = 1/2 (Small Conservative + Big Conservative)
#         - 1/2 (Small Aggressive + Big Aggressive). 	 
#
#   Rm-Rf, the excess return on the market, value-weight 
#   return of all CRSP firms incorporated in the US and 
#   listed on the NYSE, AMEX, or NASDAQ that have a CRSP
#   share code of 10 or 11 at the beginning of month t,
#   good shares and price data at the beginning of t, 
#   and good return data for t minus the one-month 
#   Treasury bill rate (from Ibbotson Associates).
#
#   See Fama and French, 1993, "Common Risk Factors 
#   in the Returns on Stocks and Bonds," Journal of 
#   Financial Economics, and Fama and French, 2014, 
#   "A Five-Factor Asset Pricing Model" for a complete
#   description of the factor returns.
#
#   Rm-Rf includes all NYSE, AMEX, and NASDAQ firms. 
#   SMB, HML, RMW, and CMA for July of year t to June 
#   of t+1 include all NYSE, AMEX, and NASDAQ stocks 
#   for which we have market equity data for December 
#   of t-1 and June of t, (positive) book equity data 
#   for t-1 (for SMB, HML, and RMW), non-missing 
#   revenues and at least one of the following: cost 
#   of goods sold, selling, general and administrative
#   expenses, or interest expense for t-1 (for SMB and
#   RMW), and total assets data for t-2 and t-1 (for 
#   SMB and CMA).
#-----------------------------------------------------
fama_french_url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_5_Factors_2x3.zip"
temp <- tempfile()
download.file(fama_french_url, temp)
fama_french <- read.table(unz(temp, 
                  "F-F_Research_Data_5_Factors_2x3.txt"),
                  skip = 3, nrows = 621) 
unlink(temp)
rm(fama_french_url, temp)
fama_french <- data.frame(month = rownames(fama_french),
                          fama_french)
rownames(fama_french) <- NULL


#-----------------------------------------------------
#   In their 2015 paper "A five-factor asset pricing model"
#   Fama and French examine the ability of various factor
#   models to explain monthly excess returns for the
#   following portfolios:
#
#       - 25 Size-B/M portfolios
#       - 25 Size-OP portfolios
#       - 25 Size-Inv portfolios
#       - 32 Size-B/M-OP portfolios 
#       - 32 Size-B/M-Inv portfolios
#       - 32 Size-OP-Inv portfolios
#
#   We should eventually take a look at these, but for
#   the moment, we'll keep things simple and look at 
#   portfolios formed on size and industry. These keeps
#   the number of parameters we need to estimate from
#   exploding.
#-----------------------------------------------------


#-----------------------------------------------------
# Portfolios Formed on Size
# Source: Ken French's Website
#-----------------------------------------------------
# DOCUMENTATION:
# 
#  Monthly Returns: July 1926 - March 2015
#
#  Portfolios: ME < 0 (not used); bottom 30%, 
#              middle 40%, top 30%; quintiles; 
#              deciles.
#
#  Construction: The portfolios are constructed at 
#                the end of each June using the 
#                June market equity and NYSE 
#                breakpoints.
#
#  Stocks: The portfolios for July of year t to June 
#          of t+1 include all NYSE, AMEX, and NASDAQ 
#          stocks for which we have market equity data 
#          for June of t.
#-----------------------------------------------------
# NOTE: These portfolios are available in both 
#       equally-weighted and value-weighted 
#       versions.
#-----------------------------------------------------
size_portfolios_url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Portfolios_Formed_on_ME.zip"
temp <- tempfile()
download.file(size_portfolios_url, temp)
size_portfolios_value <- read.table(unz(temp, 
                  "Portfolios_Formed_on_ME.txt"),
                  skip = 13, nrows = 1065, 
                  na.strings = c("-99.99", "-999")) 

size_portfolios_equal <- read.table(unz(temp, 
                  "Portfolios_Formed_on_ME.txt"),
                  skip = 1082, nrows = 1065, 
                  na.strings = c("-99.99", "-999")) 
unlink(temp)
size_names <- c("month", "Lo10", "Dec2", "Dec3", 
                "Dec4", "Dec5", "Dec6", "Dec7", 
                "Dec8", "Dec9", "Hi10")
size_portfolios_value <- size_portfolios_value[,-c(2:10)]
names(size_portfolios_value) <- size_names
size_portfolios_equal <- size_portfolios_equal[,-c(2:10)]
names(size_portfolios_equal) <- size_names
rm(size_names, temp, size_portfolios_url)

#-----------------------------------------------------
# 10 Industry Portfolios
# Source: Ken French's Website
#-----------------------------------------------------
# DOCUMENTATION:
#
#   Monthly Returns: July 1926 - December 2014
#
#   We assign each NYSE, AMEX, and NASDAQ stock 
#   to an industry portfolio at the end of June
#   of year t based on its four-digit SIC code 
#   at that time. (We use Compustat SIC codes 
#   for the fiscal year ending in calendar year 
#   t-1. Whenever Compustat SIC codes are not 
#   available, we use CRSP SIC codes for June
#   of year t.) We then compute returns from 
#   July of t to June of t+1.
#
#   1 NoDur:  Consumer NonDurables -- Food, Tobacco, Textiles, 
#             Apparel, Leather, Toys
#   2 Durbl:  Consumer Durables -- Cars, TV's, Furniture,
#             Household Appliances
#   3 Manuf:  Manufacturing -- Machinery, Trucks, Planes, 
#             Chemicals, Off Furn, Paper, Com Printing
#   4 Enrgy:  Oil, Gas, and Coal Extraction and Products
#   5 HiTec:  Business Equipment -- Computers, Software, and 
#             Electronic Equipment
#   6 Telcm:  Telephone and Television Transmission
#   7 Shops:  Wholesale, Retail, and Some Services (Laundries, 
#             Repair Shops)
#   8 Hlth:   Healthcare, Medical Equipment, and Drugs
#   9 Utils:  Utilities
#  10 Other:  Other -- Mines, Constr, BldMt, Trans, Hotels, 
#             Bus Serv, Entertainment, Finance
#-----------------------------------------------------
# NOTE: These portfolios are available in both 
#       equally-weighted and value-weighted 
#       versions.
#-----------------------------------------------------
industry_portfolios_url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios.zip"
temp <- tempfile()
download.file(industry_portfolios_url, temp)
industry_portfolios_value <- read.table(unz(temp, 
                  "10_Industry_Portfolios.txt"),
                  skip = 12, nrows = 1062, 
                  na.strings = c("-99.99", "-999")) 

industry_portfolios_equal <- read.table(unz(temp, 
                  "10_Industry_Portfolios.txt"),
                  skip = 1078, nrows = 1062, 
                  na.strings = c("-99.99", "-999")) 
unlink(temp)
industry_names <- c("month", "NoDur", "Durbl", 
                    "Manuf", "Enrgy", "HiTec",
                    "Telcm", "Shops", "Hlth", 
                    "Utils", "Other")
names(industry_portfolios_equal) <- industry_names
names(industry_portfolios_value) <- industry_names
rm(industry_names, temp, industry_portfolios_url)

#-----------------------------------------------------
# Merge the Data
#-----------------------------------------------------
all_factors <- merge(fama_french, momentum, by = "month")
all_factors <- merge(all_factors, liquidity, by = "month")

portfolios_value <- merge(industry_portfolios_value, size_portfolios_value,
                          by = "month")
portfolios_equal <- merge(industry_portfolios_equal, size_portfolios_equal,
                          by = "month")

rm(fama_french, industry_portfolios_equal, industry_portfolios_value,
   size_portfolios_equal, size_portfolios_value, momentum, liquidity)

out_equal <- merge(all_factors, portfolios_equal, by = "month") 
out_value <- merge(all_factors, portfolios_value, by = "month") 

rm(all_factors, portfolios_equal, portfolios_value)


#-----------------------------------------------------
# Output Merge Data
#-----------------------------------------------------
setwd("~/factor-choice/")
write.csv(out_value, "data_value.csv", row.names = FALSE)
write.csv(out_equal, "data_equal.csv", row.names = FALSE)

rm(out_equal, out_value)


