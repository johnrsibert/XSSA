require("XML")
first<-function()
{
doc<-xmlTreeParse("simple.xml")
doc$doc$children$FIT
doc$doc$children$FIT$nobs
doc$doc$children$FIT-> F
getNodeSet(F,"//nobs/value")->val
val
is.list(val)
length(val)
unlist(val)->unlist.val
length(unlist.val)
unlist.val[3]
as.numeric(unlist.val[3])
}

better<-functgion()
{
nodelist<-xmlToList("simple.xml")
nodelist
nodelist$Fit$nobs
nodelist$nobs
nodelist$nobs$value
}
