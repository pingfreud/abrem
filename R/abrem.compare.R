# R package 'abrem'
# Abernethy Reliability Methods
# Implementations of lifetime data analysis methods described in
# 'The New Weibull Handbook, Fifth edition' by Dr. Robert B. Abernethy.
# August 2014, Jurgen Symynck
# Copyright 2014, Jurgen Symynck
#
# For more info, visit http://www.openreliability.org/
#
# For the latest version of this file, check the Subversion repository at
# http://r-forge.r-project.org/projects/abernethy/
#
# Disclaimer:
#    The author is not affiliated with Dr. Abernethy or Wes Fulton - CEO of
#    Fulton Findings(TM) and author of the software package SuperSMITH
#-------------------------------------------------------------------------------
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# +-----------------------------------+
# |  execute this software with R:    |
# |  http://www.r-project.org/        |
# +-----------------------------------+
#
# Implementation of Life Comparison with methods suggested by Gerald G. Brown and Hebert C. Rutemiller
# as described in http://reliawiki.org/index.php/Comparing_Life_Data_Sets
# Takes beta and eta parameters for 2 distributions and returns probability P(t2 >= t1)

abrem.compare <- function(beta1,eta1,beta2,eta2){
  # Weibull 2-parameter probability density function
  weib.pdf <- function (x){
    (beta1/eta1)*((x/eta1)^(beta1-1))*(exp(1)^-((x/eta1)^beta1))
  }
  # Weibull Reliability function (2-parameter Cumulative Density Function)
  weib.cdf <- function (x){
    1 - exp(1)^-(x/eta2)^beta2 
  }
  # Composite function and integration (0 to +inf) to evaluate P (t2 >= t1)
  weib.mult <- function (x){
    weib.pdf(x)*weib.cdf(x)
  }
  dif.value <- integrate(f=weib.mult,lower = 0, upper = Inf)
  return(dif.value$value) # 
}
