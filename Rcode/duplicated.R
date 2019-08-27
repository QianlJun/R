test<- data.frame(x = c(1:10, 1:3), y=1:13)
test
test2<-test[!duplicated(test$x), ]
test2