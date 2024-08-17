rm(list=ls())

library(dplyr)
library(tidyverse)

# this needs only coverage as input!
#use name to construct the output!

#Load csv table
coverage_table <- read.csv2("A2.HPV16.coverage.csv", header= TRUE, sep= "\t")


# Transofme all needed column as numeric (might not be useful on the server)
coverage_table<-transform(coverage_table, coverage= as.numeric(as.character(coverage)), A= as.numeric(as.character(A)), 
                          G= as.numeric(as.character(G)), C= as.numeric(as.character(C)), T= as.numeric(as.character(T)), qA= as.numeric(as.character(qA)),
                          qG= as.numeric(as.character(qG)), qC= as.numeric(as.character(qC)), qT= as.numeric(as.character(qT)))


coverage_table <- coverage_table[order(coverage_table$position),]



#split df as a list of dataframes
dflist<-split(coverage_table, coverage_table$sample_id)



# This is the function that will allow to add each row for specific columns.
combine_rows_coverage <- function(data, row1, row2, row3, row4, col1, col2, col3, col4) {
  data[c(row3:row4), c(col3:col4)] <- (data[c(row1:row2), c((col1+1):col2)] * data[c(row1:row2), c(col3:col4)]+ 
                                         data[c(row3:row4), c((col1+1):col2)] * data[c(row3:row4), c(col3:col4)])
  data[c(row3:row4), c(col1:col2)] <- (data[c(row1:row2), c(col1:col2)] + data[c(row3:row4), c(col1:col2)])
  data[c(row3:row4), c(col3:col4)] <- (data[c(row3:row4), c(col3:col4)] / (data[c(row3:row4), c((col1+1):col2)]))
  data[ ,c(col3:col4)][is.na(data[ , c(col3:col4)])] <- 0
  data[-c(row1:row2), ]
}


# Move the 1000 last rows to the rows 1001 to 2000
dflist2<-lapply(dflist, function(x) {
  Start_line_top2<- 1001
  End_line_top2<- 2000
  Start_line_bottom2<- nrow(x)-999
  End_line_bottom2<- nrow(x)
  combine_rows_coverage(x, Start_line_bottom2, End_line_bottom2, Start_line_top2, End_line_top2, 4, 8, 11, 14)
  } )

# Move rows 1 to 1000 to the 1000 last rows
dflist3<-lapply(dflist2, function(x) {
  Start_line_top3<- 1
  End_line_top3<- 1000
  Start_line_bottom3<- nrow(x)-999
  End_line_bottom3<- nrow(x)
  combine_rows_coverage(x, Start_line_top3, End_line_top3, Start_line_bottom3, End_line_bottom3, 4, 8, 11, 14)
  } )




# Bind all the dataframes together, had to use this because plyr is not installed
total <- as.data.frame(do.call(rbind, dflist3), stringsAsFactors = FALSE)

#Subtract 1000 from positions to get correct positions
total <- total %>% mutate(position = position - 1000)

# Export table
write.table(total, "A2_coverage_no_overhangs.csv", sep="\t", row.names = FALSE)


