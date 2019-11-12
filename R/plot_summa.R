#' Plot predicted versus actual performance
#'
#' This is a convenience function to plot actual vs predicted performances.
#' Running this function only makes sense if the user knows the true labels and
#' have already run the calculate_performance() functions using the actual labels.
#' This function creates two plots the first plot shows the actual vs summa estimated
#' performances for each method as well as the performance of SUMMA and woc ensemble when
#' all methods are used. The second plot shows how the performance is changing if you add new
#' methods in the ensemble. In this plot for SUMMA we use the estimated performance to rank
#' methods and WOC we randomly sample methods (to make this robust we perform randomization 50 times)
#' and plot the mean performance and standart error of the mean.
#' @param summa A summa object created using the calculate_performance() function
#' @return A plot that compares actual v.s. predicted performance as well as
#' the output of the majority vote and summa ensemble
#' @examples
#' data=create_predictions(20,30,0.3,"binary")
#' summa=summa(data$predictions,"binary")
#' summa=calculate_performance(summa,data$actual_labels)
#' summa_plot(summa)
#' @export

summa_plot <- function(summa)
{

  if(summa@type=="binary")
  {
  xlabel="Estimated BA"
  ylabel="Actual BA"
  xlabel1="Number of methods"
  ylabel1="Balanced Accuracy"
  }

  if(summa@type=="rank")
  {
    xlabel="Estimated AUC"
    ylabel="Actual AUC"
    xlabel1="Number of methods"
    ylabel1="AUC"
  }
  p=ggplot2::qplot(summa@estimated_performance,summa@actual_performance)+
  ggplot2::geom_point(ggplot2::aes(colour="Individual"),size=3)+
  ggplot2::geom_abline(slope=1,intercept=0,color="red",size=1)+
    ggplot2::geom_hline(ggplot2::aes(yintercept=summa@summa,colour="SUMMA"),size=1,linetype=2)+
    ggplot2::geom_hline(ggplot2::aes(yintercept=summa@woc,colour="WOC"),size=1,
               linetype="dotdash")+
    ggplot2::scale_colour_manual("METHODS",values = c("black","red", "blue"))+
    ggplot2::scale_linetype_manual(name="METHODS",
                          values=c("solid", 2,"dotdash"))+
    ggplot2::xlab(xlabel)+
    ggplot2::ylab(ylabel)+
    ggplot2::theme(axis.text=ggplot2::element_text(size=18),
                   axis.title=ggplot2::element_text(size=20,face="bold"),
                   legend.text=ggplot2::element_text(size=18))

  plot(p)

  performance_all_estimated=summa@actual_performance
  idx=order(abs(performance_all_estimated),decreasing = TRUE)

  woc_data=data.frame(summa@woc_combined_performance,stringsAsFactors = FALSE)
 p= ggplot2::qplot(1:length(idx),summa@actual_performance[idx])+
    ggplot2::geom_point(ggplot2::aes(colour="Individual"),size=1)+
    ggplot2::geom_line(ggplot2::aes(1:length(idx),summa@summa_combined_performance,colour="SUMMA"),size=1)+
    ggplot2::stat_summary(inherit.aes = FALSE,data=woc_data,fun.data=get_SE, geom="errorbar", ggplot2::aes(n_methods,performance,colour="WOC"),size=1)+
    ggplot2::stat_summary(inherit.aes = FALSE,data=woc_data,fun.y = mean,geom = "line",ggplot2::aes(n_methods,performance),group=1,color="blue",size=1)+
    ggplot2::stat_summary(inherit.aes = FALSE,data=woc_data,fun.y = mean,geom = "point",ggplot2::aes(n_methods,performance),group=1,color="blue")+
    ggplot2::scale_colour_manual("METHODS",values = c("black","red", "blue"))+
    #ggplot2::scale_linetype_manual(name="METHODS",
    #                               values=c("solid", 2,"dotdash"))+
    ggplot2::xlab(xlabel1)+
    ggplot2::ylab(ylabel1)+
    ggplot2::theme(axis.text=ggplot2::element_text(size=18),
         axis.title=ggplot2::element_text(size=20,face="bold"),
         legend.text=ggplot2::element_text(size=18))
plot(p)


}
