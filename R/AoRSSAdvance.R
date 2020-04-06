#' @title Analysis of residual sum of squares (AoRSS) advanced
#'
#' @description Returns the test statistic and p-value corresponding to a method of comparing a number of curves called the analysis of residual sum of squares (AoRSS). The null hypothesis of test is that the K curves can be assumed coincident.
#'
#' @param data a list containing the data (in data.frame) for fitting each the individual curves. Length-Age data if type="VBG" or Weight-Length data if type="WL".The Logistic curve can be used in many cases of study, hence the data.frame can be in the first column the response variable and in the second one the points where the logistic function must be computed.
#' @param type the selected type of curve: Von Bertalanfly growth curve (type="VBG), length-weight curve (type="WL") and logistic curve (type="Logistic").
#' @param parameters the corresponding parameters according the selected type. If type="VBG", L_inf=Asymptotic average maximum body size, k=Growth rate coefficient that determines how quickly the maximum is attained, and t0=Hypothetical age at which the species has zero length. If type="WL", a=Allometric growth parameter and b=Scaling constant. If type="Logistic", x50=x-value of the sigmoid's midpoint of the logistic function and xd=Minus the inverse of the logistic growth rate (steepness of the curve).
#' @details See details of AoRSS.test() function. The advance of the current function is that the residuals are computed inside the function and the user only needs to provide de parameters of the individual and pooled curves. Of course, this function implement a limit type of curves: Von Bertalanfly growth curve and length-weight curve.\itemize{
#' \item Length-Weight relationship
#'
#' \deqn{W=a*L^b}
#'
#' where L is the vector of lengths by age, a is the allometric growth parameter, b scaling constant, and W is the age weight vector.
#' \item Von Bertalanffy Growth Model
#'
#'\deqn{L(x)=L_inf*(1-exp(-k*(x-t0)))}
#'
#' where L_inf is the asymptotic average maximum body size, t0 is hypothetical age at which the species has zero length, and k is growth rate coefficient.
#'
#' \item Logistic function \deqn{L(x)=1/(1+exp((x-x50)/xd))} where x50 is the x-value of the sigmoid's midpoint of the logistic function, and xd is minus the inverse of the logistic growth rate (steepness of the curve).}
#' @return The value of the F-statistic and the corresponding p-value.
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' }
#' @examples
#' # An example based on the age length data relating
#' # to Pacific hake with separate data for both males
#' # and females (Example Table 9.3 of Haddon 2011,
#' # Quantitative Methods of Fisheries). The question
#' # is whether the male and female Pacific hake exhibit
#' # different growth throughout their lives. This is, we
#' # testing if the growth curves for males and females
#' # are coincident or not?
#' # Now, we create the argument parameters.
#' # The values for each curve must be in the same order
#' # as the data in our list, and the last value will
#' # be for the pooled curve.
#' # In our case the first value corresponds to Females,
#' # the second to Males, and the last one is for the
#' # pooled curve.
#' Linf=c(61.2332329172338, 55.9779612903168, 59.2937504702457)
#' t0=c(-0.057255722705518,0.171338088033198,0.010431548715808)
#' k=c(0.29625474003956, 0.385585512821354,0.32048121909873)
#' parameters=list(Linf=Linf,t0=t0,k=k)
#' # Note that in the data.frames the age and the length must be
#' # in columns with names "Age" and "Lenght".
#' Female=data.frame(Age=c(1.0,2.0,3.3,4.3,5.3,6.3,7.3,8.3,9.3,
#' 10.3,11.3,12.3,13.3),
#' Length=c(15.40,28.03,41.18,46.20,48.23,50.26,51.82,54.27,56.98,
#' 58.93,59.00,60.91,61.83))
#' Male=data.frame(Age=c(1.0,2.0,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3),
#' Length=c(15.40,26.93,42.23,44.59,47.63,49.67,50.87,52.30,54.77,
#' 56.43,55.88))
#' data=list()
#' data[[1]]=Female
#' data[[2]]=Male
#' AoRSS.test.advanced(data,parameters,type="VBG")
#'
#' # For type="WL" is the same, but the parameters
#' # are a and b, and the data.frame must have Length
#' # and Weight (the colnames of each data.frame must be
#' # exactly such names).
#'
#' # For type="Logistic". We have two curves of maturity that we nees to compare.
#' x50=c(3,3.1,2.95)
#' xd=c(-0.5,-0.62,0.51)
#' parameters=list(x50=x50,xd=xd)
#' data=list()
#' ages=1:5; Mat1=c(0.008162571, 0.057324176, 0.310025519, 0.768524783, 0.960834277)
#' Mat2=c(0.01212843,0.08317270,0.40131234,0.83201839,0.97340301)
#' data[[1]]=data.frame(Mat1,ages)
#' data[[2]]=data.frame(Mat2,ages)
#' AoRSS.test.advanced(data,parameters,type="Logistic")
#' @export

AoRSS.test.advanced<-function(data,parameters,type){
  if(type=="VBG"){

  K=length(parameters$Linf)-1

  RSSi=1:K
  Linf=parameters$Linf
  t0=parameters$t0
  k=parameters$k

  for (i in 1:K){
  datai=data[[i]]
  est=Linf[i]*(1-exp(-k[i]*(as.numeric(datai$Age)-t0[i])))

  RSSi[i]=sum((as.numeric(datai$Length)-est)^2)
  }

  # TOTAL
  age=0;lengthd=0
  for(i in 1:K){
    datai=data[[i]]
    age=c(age,datai$Age)
    lengthd=c(lengthd,datai$Length)
  }
  lengthd=lengthd[-1]
  age=age[-1]

  estT=Linf[K+1]*(1-exp(-k[K+1]*((as.numeric(age))-t0[K+1])))

  RSSp=sum((as.numeric(lengthd)-estT)^2)

  N=length(lengthd)
  return(AoRSS.test(RSSi,RSSp,K,N))
  }

  if(type=="WL"){

    K=length(parameters$a)-1

    RSSi=1:K
    a=parameters$a
    b=parameters$b


    for (i in 1:K){
      datai=data[[i]]
      est=a[i]*(as.numeric(datai$Length))^b[i]

      RSSi[i]=sum((as.numeric(datai$Weight)-est)^2)
    }

    # TOTAL
    weightd=0;lengthd=0
    for(i in 1:K){
      datai=data[[i]]
      weightd=c(weightd,datai$Weight)
      lengthd=c(lengthd,datai$Length)
    }
    lengthd=lengthd[-1]
    weightd=weightd[-1]

    estT=a[K+1]*(as.numeric(lengthd))^b[K+1]

    RSSp=sum((as.numeric(weightd)-estT)^2)

    N=length(lengthd)
    return(AoRSS.test(RSSi,RSSp,K,N))
  }

  if(type=="Logistic"){
    Logistic<-function(x,x50,xd){
      if(x[1]==0){x[1]<-0.0001}
      Lo<-1/(1+exp((x-x50)/xd))
      return(Lo)
    }

    K=length(parameters$x50)-1

    RSSi=1:K
    x50=parameters$x50
    xd=parameters$xd


    for (i in 1:K){
      datai=data[[i]]
      est=Logistic(datai[,2],x50[i],xd[i])

      RSSi[i]=sum((as.numeric(datai[,1])-est)^2)
    }

    # TOTAL
    dat1=0;dat2=0
    for(i in 1:K){
      datai=data[[i]]
      dat1=c(dat1,datai[,1])
      dat2=c(dat2,datai[,2])}
    dat1=dat1[-1]
    dat2=dat2[-1]

    estT=Logistic(dat2,x50[K+1],xd[K+1])

    RSSp=sum((as.numeric(dat1)-estT)^2)

    N=length(dat1)
    return(AoRSS.test(RSSi,RSSp,K,N))
  }

}
