#读取数据集
dataset <- read.table("D:/code/r_language_code/learn_task/t_text.txt",header =T)
dataset <- data.frame(dataset)

#单样本检验封装函数
one_simple_t_test <- function(test_var1,mu,t_dataset){
  result <- t.test(t_dataset[,test_var1], alternative = "two.side", mu = mu )

  attr(result$statistic,"names") <- NULL  #检验结果t值
  t_value <- result$statistic

  attr(result$parameter,"names") <- NULL  #检验结果df
  df <- result$parameter

  p_value <- result$p.value  #检验结果p值
  result_list <- list(t_value= t_value , df =df , p_value = p_value )
  return(result_list)
}

one_simple_t_test("x",14.02,dataset)


#读取数据集
library(readxl)
dataset2  <- read_excel("D:/code/r_language_code/learn_task/paird_t.xlsx" )
dataset2 <- data.frame(dataset2)
t.test(dataset2$x1,dataset2$x2,paired=TRUE)

#配对样本函数封装
paired_t_test <- function(test_var1,test_var2,dataset){
  result <- t.test(dataset[,test_var1],dataset[,test_var2],paired=TRUE)

  attr(result$statistic,"names") <- NULL  #检验结果t值
  t_value <- result$statistic

  attr(result$parameter,"names") <- NULL  #检验结果df
  df <- result$parameter

  p_value <- result$p.value  #检验结果p值

  conf_int_upper <- result$conf.int[1]  #检验结果置信区间
  conf_int_lower <- result$conf.int[2]

  attr( result$estimate,"names") <-NULL  #检验结果均值
  estimate <- result$estimate
  result_list <- list(t_value= t_value , df =df , p_value = p_value, conf_int_upper=conf_int_upper,
                      conf_int_lower = conf_int_lower,estimate=estimate)
  
}
result <- paired_t_test("x1","x2",dataset2)

