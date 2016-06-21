

# library(pander)


#######################################
## Crosstab for basic count statitics
#######################################
#' Creates a Crosstable with basic counts
#'
#'
#'
#' @param row vector containing factors for row
#' @param col vector containing factors for columns
#' @param  margin index number (1 for rows, etc.)
#' @param  dec decimal points (default =1)
#' @param total sum of columns and rows (default= FALSE)
#'
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatientCohort<-c("CohortA","CohortA","CohortA","CohortB","CohortB","CohortB","CohortB","CohortB","CohortB","CohortC","CohortC","CohortC","CohortC")
#' PatientClass<-c("ClassA","ClassB","ClassB","ClassC","ClassA","ClassC", "ClassA","ClassB","ClassA","ClassC","ClassA","ClassA","ClassA")
#' t<-ctab(row=PatientCohort,col=PatientClass,margin=1,total=T,perc=T)
#'
#'@export
ctab <- function(row, col, margin = 1, dec = 1, percs = FALSE, total = FALSE){
  #tab <- as.table(table(row[[1]],col[,1])) #stb 20160607: row war data.frame jetzt vector
  if (is.data.frame(row)==T){row<-row[[1]]}
  if (is.data.frame(col)==T){col<-col[[1]]}
  tab <- as.table(table(row,col)) #stb 20160607: row war data.frame jetzt vector
  innermargin<-1
  if (margin==2){innermargin<-2}
  ptab <- signif(prop.table(tab, margin = innermargin)*100, dec)

  if (percs){

    z <- matrix(NA, nrow = nrow(tab), ncol = ncol(tab), byrow = TRUE)
    for (i in 1:ncol(tab)) z[,i] <- paste(tab[,i], ptab[,i], sep = " ")
    rownames(z) <- rownames(tab)
    colnames(z) <- colnames(tab)

    if (margin == 1 & total){
      rowTot <- paste(apply(tab, 1, sum)," (", apply(ptab, 1, sum),")", sep = "")
      z <- cbind(z, Gesamt = rowTot)
    } else if (margin == 2 & total) {
      colTot <- paste(apply(tab, 2, sum)," (", apply(ptab, 2, sum),")", sep = "")
      z <- rbind(z,Gesamt = colTot)
    }else if (margin == 3 & total) {
      tab_t<-rbind(tab,Gesamt=apply(tab, 2, sum))
      ptab_t <- signif(prop.table(tab_t, margin = 1)*100, dec)

      z <- matrix(NA, nrow = nrow(tab_t), ncol = ncol(tab_t), byrow = TRUE)
      for (i in 1:ncol(tab_t)) z[,i] <- paste(tab_t[,i], " (",ptab_t[,i],")", sep = "")
      rownames(z) <- rownames(tab_t)
      colnames(z) <- colnames(tab_t)
      rowTot <- paste(apply(tab_t, 1, sum)," (", signif(apply(ptab_t, 1, sum),1),")", sep = "")
      z <- cbind(z,Gesamt = rowTot)
    }
  }
  else {
    z <- table(row, col)

  }
  return(z)

}
#############################
#
# Crosstab with 2 dimensions an one variable for unique counting
#
#############################
#' Creates a Crosstable in Markdown format with 2 dimensions and one variable for a unique counting
#'
#'
#'
#' @param row vector containing factors for row
#' @param col vector containing factors for columns
#' @param data the dataframe to get the values from
#' @param margin index number (1 for rows, etc.)
#' @param dec decimal points (default =1)
#' @param total sum of columns and rows (default= TRUE)
#' @param percs if TRUE percentage of row or column are displayed (default=TRUE)
#' @param split.tables where to split wide tables to separate tables. The default value is 300 characters
#' @param margin index, or vector of indices to generate margin for
#' @param table.continues: string (default: 'Tabelle wird fortgesetzt') passed to pandoc.table to be used as caption for long (split) without a use defined caption
#' @param sortcounts if TRUE rows will be sorted descendent by counts
#' @param split.cells where to split cells' text with line breaks. Default to 30, to disable set to Inf. Can be also supplied as a vector, for each cell separately (if length(split.cells) == number of columns + 1, then first value in split.cells if for row names, and others are for columns). Supports relative (percentage) parameters in combination with split.tables.
#' @param ucol.name title of counted values default is "Patienten [n/(%)]"
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatientCohort<-c("CohortA","CohortA","CohortA","CohortB","CohortB","CohortB","CohortB","CohortB","CohortB","CohortC","CohortC","CohortC","CohortC")
#' PatientClass<-c("ClassA","ClassB","ClassB","ClassC","ClassA","ClassC", "ClassA","ClassB","ClassA","ClassC","ClassA","ClassA","ClassA")
#' df<-as.data.frame(cbind(PatientCohort,PatientClass))
#' tab<-bbcounttab(row="PatientCohort",col="PatientClass",data=df,total=T,perc=T,dec=2,ucol.name="Classes")
#' tab
#'
#'@export
bbcounttab<-function(row=NULL,col=NULL,data,uniquecounts=NULL,caption=NULL,percs=TRUE,total=TRUE,dec=1,split.tables=300,margin=3,
                     table.continues="Tabelle wird fortgesetzt" ,sortcounts=TRUE,split.cells=30,ucol.name="Patienten [n/(%)]"){
  pander::panderOptions('table.continues',table.continues)
  pander::panderOptions('table.caption.prefix','Table: ')
  # panderOptions('table.caption.prefix','Tabelle: ')
  # In case we have only one dimension for our table (simplifies things)
  data<-bbsurvr::drop.levels(data)
  if ((is.null(row)==TRUE) | (is.null(col)==TRUE)){
    fac<-coalesce(row,col)
    if(is.null(uniquecounts)==TRUE){daten<-data[,c(fac)]}else{daten<-unique(data[,c(fac,uniquecounts)])}

    daten<-cbind(daten,freq=1)
    if(sortcounts==TRUE){
      stab<-sort(xtabs(freq~daten[,fac],daten),decreasing = TRUE)
      ptab<-sort(round(prop.table(stab)*100,dec),decreasing = TRUE)
    }else{
      stab<-xtabs(freq~daten[,fac],daten)
      ptab<-round(prop.table(stab)*100,dec)
    }
    mlen<-coalesce(nrow(stab),length(names(stab)))
    z <- matrix(NA, ncol = mlen,  byrow = TRUE)
    colnames(z) <- coalesce(rownames(stab),names(stab))
    if (is.null(dim(stab))){
      z[1,] <- paste(stab," (", signif(ptab,3),")", sep = "")
    }else{
      z[1,] <- paste(apply(stab, 1, sum)," (", signif(apply(ptab, 1, sum),3),")", sep = "")}
    z<-cbind(z,Gesamt = paste(sum(stab)," (", sum(ptab),")", sep = ""))
    if(is.null(row)){
      tab<-t(z)
      names(tab)<-rownames(stab)
      if(is.null(ucol.name)==FALSE){colnames(tab)<-ucol.name}
    }else{
      tab<-z}
  }
  else{
    if(is.null(uniquecounts)==TRUE){daten<-data[,c(row,col)]}else{daten<-unique(data[,c(row,col,uniquecounts)])}

    tab<-ctab(row=daten[,row],col=daten[,col],percs=percs,total=total,dec=dec,margin=margin)
    # alignment der tabelle
    #set.alignment(default = "centre", row.names = "left")
  }
  set.alignment(default = "centre", row.names = "left")
  pander::pandoc.table(tab,split.tables=split.tables,split.cells=split.cells,style="multiline",caption=caption)
}

#############################
#
# Crosstab with 1 dimensions and summary of one variable (Mean,Sd,median 1,2nd quantile, Counts)
#
#############################

#' Creates a markdown Crosstable with 1 dimensions and summary of one variable (Mean,Sd,median 1,2nd quantile, Counts)
#'
#' @import dplyr
#'
#' @param data the dataframe to use
#' @param var the variable to count
#' @param group to group "var" by
#' @param caption of the table
#' @param total sum of columns and rows (default= TRUE)
#' @param dec decimal points (default =1)
#' @param split.tables where to split wide tables to separate tables. The default value is 300 characters
#' @param margin index number (1 for rows, etc.)
#' @param table.continues: string (default: 'Tabelle wird fortgesetzt') passed to pandoc.table to be used as caption for long (split) without a use defined caption
#' @param groupname column header for groups
#' @param orderby the column to order by

#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatientWeight<-c(66.5,82.2,74.8,70.2,95.7,55.8,59.2,77.2,75.0,75.0,102.8,62.8,65.3)
#' PatientClass<-c("ClassA","ClassB","ClassB","ClassC","ClassA","ClassC", "ClassA","ClassB","ClassA","ClassC","ClassA","ClassA","ClassA")
#' df<-as.data.frame(cbind(PatientWeight,PatientClass))
#' tab<-bbsummarizetab(data=df,var="PatientWeight",group="PatientClass",caption="Summary of patients weight",total=T,dec=2,groupname="Classes")
#' tab
#'
#'@export
bbsummarizetab<-function(data,var,group=NULL,caption=NULL,total=TRUE,dec=1,split.tables=300,margin=3,
                         table.continues="Tabelle wird fortgesetzt",groupname=NULL,orderby=NULL ){
 # require(dplyr)# hier mus as.numerich verwenbdet werden wg bug: https://github.com/hadley/dplyr/issues/893

  data<-subset(data,is.na(data[,var])==FALSE)

  colnames(data)[colnames(data)==var]<-"uvar"
  data$uvar<-as.numeric(data$uvar)
  tab0<-data%>%summarise( Min=min(uvar)
                          ,Q25= quantile(uvar,probs=0.25)
                          ,Mittel = mean(uvar)
                          ,SD = sd(uvar)
                          ,Median = median(uvar)
                          ,Q75 = quantile(uvar,probs=0.75)
                          ,Max=max(uvar)
                          ,"Anzahl\nWerte"=length(uvar))


  if (is.null(group)==FALSE){

    colnames(data)[colnames(data)==group]<-"ugroup"

    if(!is.null(orderby)){
      colnames(data)[colnames(data)==orderby]<-"ordervar"
      data<-data%>%arrange(ordervar)
      tab<- data%>%group_by("group"=ugroup,ordervar)%>%summarise( Min=min(uvar)
                                                                  ,Q25= quantile(uvar,probs=0.25)
                                                                  ,Mittel = mean(uvar)
                                                                  ,SD = sd(uvar)
                                                                  ,Median = median(uvar)
                                                                  ,Q75 = quantile(uvar,probs=0.75)
                                                                  ,Max=max(uvar)
                                                                  ,Anz=length(uvar))%>%ungroup()%>%arrange(ordervar)
      tab<-tab%>%select(group,Min,Q25,Mittel,SD,Median,Q75,Max,"Anzahl\nWerte"=Anz)
      #dataorder<-tab%>%select(group)%>%distinct()
      #levels(tab$group)<-dataorder[,1]
    }else{

      tab<- data%>%group_by("group"=ugroup)%>%summarise( Min=min(uvar)
                                                         ,Q25= quantile(uvar,probs=0.25)
                                                         ,Mittel = mean(uvar)
                                                         ,SD = sd(uvar)
                                                         ,Median = median(uvar)
                                                         ,Q75 = quantile(uvar,probs=0.75)
                                                         ,Max=max(uvar)
                                                         ,"Anzahl\nWerte"=length(uvar))
    }



    if (total==TRUE){
      #names(tab0)<-names(tab)
      tab0<-cbind(group,tab0)
      levels(tab0$group)<-c(levels(tab0$group),"Gesamt")
      tab0[1,1]<-"Gesamt"
      tab<-rbind(tab,tab0)
    }
    if(is.null(groupname)==FALSE){colnames(tab)[1]<-groupname }else{colnames(tab)[1]<-"Gruppe"}
  }
  else {tab<-tab0[,1:length(tab0)]}
  pander::panderOptions('table.continues',table.continues)
  pander::panderOptions('table.caption.prefix','Table: ')
  if(is.null(group)==FALSE){pander::emphasize.strong.cols(1)}
  pander::set.alignment(default = "centre", row.names = "left")
  pander::pandoc.table(tab,split.tables=split.tables,style="multiline",caption=caption)

}
#########################################################
# BArplot with ggplot
#########################################################
#' Streamlines the creation of a barplot for primitive counting  using ggplot2
#'
#' @import ggplot2
#'
#' @param data the dataframe to use
#' @param factor to group the data by
#' @param uniquecounts the variable to count
#' @param countslab lable of the y-axis
#' @param xlab label of x-axis
#' @param xrotate if true labels of x-axis are rotation by 90°
#' @param facet the factor to facet the plot by
#' @param prop if true the proportions will be ploted
#' @param xaxisorder vector of values to sort the values of "factor"
#' @param horizontal if ture a horizontal barplot is ploted
#' @param stackpar if the name of a variable is given a stacked barplot will be generated
#' @param facetncol the amount of column the facet plot is generated to
#' @param stacktitle the title of the stacked variable default="Gruppe"
#' @param facetscales the scale = "free"
#' @param cex.datalabel fontsize of datalabels default=2,
#' @param datalabel if true the data will be labeld default=TRUE

#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatID<-seq(1,13)
#' PatientWeightClass<-c("61-70","41-50","61-70","41-50","61-70","71-80","81-90","81-90","71-80","81-90","81-90","91-100","61-70")
#' PatientClass<-c("ClassA","ClassB","ClassB","ClassC","ClassA","ClassC", "ClassA","ClassB","ClassA","ClassC","ClassA","ClassA","ClassA")
#' df<-as.data.frame(cbind(PatID,PatientWeightClass,PatientClass))
#' p1<-bbbarplot(data=df,"PatientWeightClass","PatID",countslab="Anzahl", xlab=NULL,xrotate=FALSE,stackpar="PatientClass",stacktitle="PatientClass",datalabel=TRUE)
#' p1
#'
#'@export
bbbarplot<-function(data, factor,uniquecounts,countslab="Anzahl", xlab=NULL,xrotate=FALSE,facet=NULL,prop=FALSE,xaxisorder=NULL,
                      horizontal=FALSE,stackpar=NULL,facetncol=2,stacktitle="Gruppe",facetscales = "free",cex.datalabel=2,
                      datalabel=TRUE
){
  #require(ggplot2)

  #require(RColorBrewer)
  hjustdat=0.5 # default horizontal adjustment for datalabels
  vjustdat=-.2
  #names(data[,factor])<-"ufact"
  if (is.null(xlab))
  {xlab<-factor}
  colnames(data)[colnames(data)==uniquecounts]<-"ucountdat"
  if(prop==FALSE){
    if(is.null(facet)){
      if(is.null(stackpar)){
        imgdat<-data%>%group_by_(xfac=factor)%>%summarise(ycount=n_distinct(ucountdat))
      }else{
        imgdat<-data%>%group_by_(xfac=factor,stack=stackpar)%>%summarise(ycount=n_distinct(ucountdat))
      }

    }else{
      if(is.null(stackpar)){
        imgdat<-data%>%group_by_(facet=facet,xfac=factor)%>%summarise(ycount=n_distinct(ucountdat))
      }else{
        imgdat<-data%>%group_by_(facet=facet,xfac=factor,stack=stackpar)%>%summarise(ycount=n_distinct(ucountdat))
      }

    }
  }else{
    if (countslab=="Anzahl"){countslab<-"Prozent"}
    if(is.null(facet)){
      if(is.null(stackpar)){
        imgdat<-data%>%group_by_(xfac=factor)%>%summarise(Anzahl=n_distinct(ucountdat))%>%mutate(ycount=round(100/length(unique(data$ucountdat))*Anzahl,1))
      }else
      {

        imgdat<-data%>%group_by_(xfac=factor,stack=stackpar)%>%summarise(Anzahl=n_distinct(ucountdat))%>%mutate(ycount=round(100/sum(Anzahl)*Anzahl,2))

      }
    }else{
      # Das ist noch nicht gelöst
    }
  }


  if(is.null(stackpar)){
    textpos<-position_dodge(width=0.9)
    bp<-ggplot(imgdat, aes(xfac, ycount))+geom_bar(stat="identity",fill = "#004992", colour = "darkgrey", alpha = 0.8,position=position_dodge())
  } else{
    nstack<-n_distinct(imgdat$stack)
    if(nstack==2){
      cols<-c("#E41A1C" ,"#4D74AB")}
    else{
      cols<-bbsurvr::getBBColors(n_distinct(imgdat$stack))}
    textpos<-position_stack()
    bp<-ggplot(imgdat, aes(xfac, ycount,fill = stack))+geom_bar(stat="identity", colour = "darkgrey", alpha = 0.8
                                                                ,position=position_stack())+scale_fill_manual(values = cols,
                                                                                                              name=stacktitle)
    #+scale_fill_brewer(palette = 12)
  }

  bp<-bp+theme_minimal(base_size = 11, base_family = "Helvetica") + theme(axis.title = element_text(vjust=0.1),axis.title.y=element_text(vjust=0.5))+labs(x=xlab,y=countslab)
  if(!is.null(facet)){bp<-bp+facet_wrap(~facet,ncol=facetncol,scales=facetscales)}
  if(xrotate){bp<-bp+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}


  if(is.null(xaxisorder)==FALSE){# sortierung des x-achse

    bp<-bp+scale_x_discrete(limits=xaxisorder)
  }
  if(horizontal==TRUE){
    bp<-bp+coord_flip()
    hjustdat=-.2
    vjustdat=.4

  }
  if(datalabel==TRUE){
    bp<-bp+geom_text(aes(y=ycount, label=round(ycount,2)), position= textpos,hjust=hjustdat, vjust=vjustdat, color="black",size=cex.datalabel)
  }
  if(!is.null(stackpar)){
    bp<-bp+theme(legend.position="bottom",legend.key.size = unit(0.3, "cm"))

  }
  return (bp)
}
##################################################################
# Legt den Pfad zu den Reporttemplates als Globale Variable an
##################################################################
#' Set the Path to Report Template Folder as global variable
#'
#'
#'
#' @param x path path to folder

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bbSetTplPath("c:/temp/Test")
#'
#'@export
bbSetTplPath<-function(x){

  TplPath<<-x
  return(TplPath)

}
##################################################################
# Gibt unter Verwendung der Globaken Variablen TplPath das Reporttemplate zurück
##################################################################
#' Get the full Path for a Report Template based on the predefined Path by "bbSetTplPath"
#'
#'
#'
#' @param x Name of .Rmdd or .Rnw template

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' report<-bbGetTpl("c:/temp/Test.Rmd")
#'
#'@export
bbGetTpl<-function(x){


  return(paste(TplPath,"/",x,sep=""))

}

##################################################################
# Wrapperfuncction for knit_expand
##################################################################
#' Wrapper for knit_expand
#'
#'
#'
#' @param file pthe path to file

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bbreport("c:/temp/Test.Rmd")
#'
#'@export
bbreport<-function(file,...){

  scr<-knitr::knit_expand(file,...)
  return(knitr::knit(text = unlist(scr), quiet = TRUE))

}


##################################################################
# Wrapperfuncction for german dates
##################################################################
#' Formats a Datestring to German Date format (dd-mm-yyyy)
#'
#'
#'
#' @param x date string

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bbdedate("2016-11-05")
#'
#'@export
bbdedate<-function(x){

  return(strftime(x,"%d.%m.%Y"))
}
##################################################################
# Wrapperfuncction for latex '\\cleardoublepage' tag withhin markdown to control pagination in pdf output
##################################################################
#' Wrapperfuncction for latex '\\cleardoublepage' tag withhin markdown to control pagination in pdf output
#'
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bbclearpage()
#'
#'@export
bbclearpage<-function(){

  return(paste("\\cleardoublepage \n"))
}


#########################################################
# histogram with ggplot
#########################################################
#' Streamlines plotting of a histogramm using ggplot2
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param data dataframe
#' @param factor the column to count
#' @param uniquecounts the column that determines the unique id
#' @param countslab  label of y axis default="Anzahl"
#' @param xlabel label of x axis default=NULL
#' @param bin bandwidth default=1
#' @param facet if given the plot wil be printed as facets default=NULL
#' @param dens if TRUE a density plot wil be performed default=FALSE

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatID<-seq(1,50)
#' PatientWeight<-runif(50, min=40, max=100)
#' PatientWeightClass<-c("61-70","41-50","61-70","41-50","61-70","71-80","81-90","81-90","71-80","81-90","81-90","91-100","61-70")
#' PatientClass<-c("ClassA","ClassB","ClassB","ClassC","ClassA","ClassC", "ClassA","ClassB","ClassA","ClassC","ClassA","ClassA","ClassA")
#' PatientClass2<-seq(1,50,50)

#' df<-as.data.frame(cbind(PatID,PatientWeight,PatientClass))
#' p1<-bbgghistplot(data=df,"PatientWeight","PatID",countslab="Anzahl",bin=10)
#' p1
#'
#'@export
bbgghistplot<-function(data, factor,uniquecounts,countslab="Anzahl", xlabel=NULL,bin=1,facet=NULL,dens=FALSE){
  #require(ggplot2)
  if (is.null(xlabel))
  {xlabel<-factor}
  if(is.null(facet)){hdaten<-unique(data[,c(factor,uniquecounts)])
  names(hdaten)<-c("hfactor","hcounts")}else{
    hdaten<-unique(data[,c(factor,uniquecounts,facet)])
    names(hdaten)<-c("hfactor","hcounts","hfacet")
  }
  hdaten$hfactor<-as.numeric(hdaten$hfactor)# hier mus as.numerich verwenbdet werden wg bug: https://github.com/hadley/dplyr/issues/893

  mediane<-if (is.null(facet)){hdaten%>%summarise(med=median(hfactor,na.rm=TRUE))
  }else{hdaten%>%group_by(hfacet)%>%summarize(med=median(hfactor,na.rm=TRUE))}


  if(!dens){
    g<-ggplot(hdaten, aes(x=hfactor)) + geom_histogram(binwidth=bin, colour="darkgrey", fill="#004992", alpha = 1/2)
  }
  else{
    g<-ggplot(hdaten, aes(x=hfactor)) + geom_histogram(aes(y=..density..),binwidth=bin, colour="darkgrey", fill="#004992", alpha = 1/2)+geom_density(colour="darkgrey",alpha=.05,fill="red")

    countslab="Dichte"
  }
  g<-g+theme_minimal(base_size = 10) +theme(axis.title = element_text(vjust=0.1))+labs(x=xlabel,y=countslab)

  g<-g+ geom_vline(data=mediane,aes(xintercept=med), color="darkred", linetype="dashed", size=1)

  if(is.null(facet)==FALSE){g<-g+facet_wrap(~hfacet,ncol=2)}

  #g<-g+ geom_vline(aes(xintercept=median(hfactor, na.rm=T)), color="darkred", linetype="dashed", size=1)
  return(g)

}

#########################################################
# boxplot with ggplot
#########################################################
#' boxplot with ggplot
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param data dataframe
#' @param factor the factor to count
#' @param group the group
#' @param uniquecounts the unique counts
#' @param ylabel the label of y axis default="Anzahl"
#' @param xlabel the label of x axis default=NULL
#' @param xrotate if TRUE the x and y axis will switch default=FALSE
#' @param facet if a facet is given plot will be faceted default=NULL

#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' PatID<-seq(1,50)
#' PatientWeight<-runif(50, min=40, max=100)
#' PatientClass<-rep(c("ClassA","ClassB","ClassC","ClassB","ClassC"),10)
#' df<-as.data.frame(cbind(PatID,PatientWeight,PatientClass))
#' df$PatientWeight<-as.double(df$PatientWeight)
#' p1<-bbggboxplot(data=df,factor="PatientWeight",group="PatientClass",uniquecounts="PatID",ylabel="Anzahl",bin=1)
#' p1
#'
#'@export
bbggboxplot<-function(data, factor,group ,uniquecounts,ylabel="Anzahl", xlabel=NULL,xrotate=FALSE,
                      facet=NULL){
  #require(ggplot2)
  #require(dplyr)
  if (is.null(xlabel))
  {xlabel<-factor}
  bdaten<- if (is.null(facet)){unique(data[,c(factor,group,uniquecounts)])} else {unique(data[,c(factor,group,facet,uniquecounts)])}
  names(bdaten)<-if (is.null(facet)){c("bfactor","bgroup","bcounts")}else{c("bfactor","bgroup","bfacet","bcounts")}
  #mediane<-if (is.null(facet)){ddply(bdaten,.(bgroup),summarize,med=median(bfactor))}else{ddply(bdaten,.(bfacet,bgroup),summarize,med=median(bfactor))}
  mediane<-if (is.null(facet)){bdaten%>%group_by(bgroup)%>%summarise(med=median(bfactor))}else{bdaten%>%group_by(bfacet,bgroup)%>%summarise(med=median(as.numeric(bfactor)))}
  g<- ggplot(bdaten, aes(x=bgroup , y=bfactor))
  if (!is.null(facet)){g<-g+facet_wrap(~bfacet)}
  g<-g+ geom_boxplot( colour="black", fill="#004992", alpha = 1/2)+theme_minimal(base_size = 10) +
    theme(axis.title = element_text(vjust=0.1))+labs(x=xlabel,y=ylabel)
  if(xrotate){g<-g+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}

  g<-g+geom_text(data = mediane, aes(y = med, label = round(med,2)),size = 3, vjust = -0.5)
  return(g)
}

#########################################################
#Gibt einen "Count Distinct" für die Spalte eines Dataframes zurück Achtung: NA wird ausgeschlossen
#########################################################
#' Performs a distinct count of Items in a vector (omits NA). Basically a wrapper for length() and unique()
#'
#' @param data a vector
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bbDistinctCount(c(12,15,12,17,19,22,24,23,23,66,78,65,34,NA))
#'
#'@export
bbDistinctCount<-function(data){

  return(length(na.omit(unique(data))))

}

#########################################################
# Prüft, ob ein Dataframe Daten hat
# Return: bool
#########################################################
#' Checks if a Dataframe is empty or not. If TRUE
#'
#' @param data a dataframe
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' df<-as.data.frame(c("a","b","c"))
#' bbChkDataFrame(df)
#'
#'@export
bbChkDataFrame<-function(data){
  retval<-FALSE
  if(length(unique(data[,1]))>0){
    retval<-TRUE
  }
  return(retval)

}

#######################################################################################################
#######################################################################################################
#' Multiple plot function
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param cols:   Number of columns in layout
#' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' @param plotlist: vector of plots
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#'
#'@export
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}

# Knitr Hook function für den Abbildunsnamen
#
# ggp
# knit_hooks$set(abbLab = function(before, options, envir) {
#   if (before){
#     if (exists("AbbildungsCounter")==FALSE){AbbildungsCounter<<-1}
#
#     AbbLab<<-paste("Abbildung_",AbbildungsCounter,sep="")
#     options$label<-AbbLab
#     AbbildungsCounter<<-AbbildungsCounter+1
#   }
# })

##########################################################################################################
#
# Plot der rgb Praxenkarte
#
##########################################################################################################
#' Plot of Map
#'
#'  #@import rgdal
#' # @import maptools
#' # @import ggplot2
#'
#' @param x dataframew with long and lat
#' @param shapefile data source name (interpretation varies by driver — for some drivers, dsn is a file name, but may also be a folder)
#' @param layer layer name (varies by driver, may be a file name without extension)
#'
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' readOGR(dsn=paste(c(GetRDirPath(),"/Data/ShapeGER/"),collapse=""), layer="DEU_adm1",verbose=FALSE)
#'
#'
#'@export
#'
oncosys_map<-function(x,count.breaks=waiver(),shapefile,layer){

  #require("rgdal") # requires sp, will use proj.4 if installed
 # require("maptools")
  #require("ggplot2")

  pats<-x

  # Plot Theme
  theme_opts <- list(theme(panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.background = element_blank(),
                           legend.key = element_rect(fill="white"),
                           legend.title = element_text(size = rel(1.5)),
                           legend.text = element_text(size = rel(1)),
                           plot.background = element_rect(fill="white"),
                           panel.border = element_blank(),
                           axis.line = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           plot.title = element_text(size=22)))


  deu = readOGR(dsn=shapefile, layer=layer,verbose=FALSE)
  deu@data$id=rownames(deu@data)
  deu_robin <- spTransform(deu, CRS("+proj=merc"))
  deu_df_robin <- fortify(deu_robin)
  deu.df = fortify(deu)

  #Praxendaten aufbereiten
  pats_robin_df <- project(cbind(pats$Long,pats$Lat ), "+proj=merc")
  pats_robin_df <- as.data.frame(pats_robin_df)
  names(pats_robin_df) <- c("LONGITUDE", "LATITUDE")
  pats_robin_df$Patienten <- pats$anz_anm

  #PLOT

  #ggplot(deu.df) +
  ggplot(deu_df_robin) +
    aes(long,lat,group=group) +
    geom_polygon(colour= "#4B4B4B", fill="#DEDEDE") +
    geom_path(color="#4D74AB") + theme_opts +
    geom_point(data=pats_robin_df, aes(LONGITUDE ,LATITUDE, group=NULL, fill=NULL, size=Patienten), color="#4D74AB", alpha=I(8/10))+
    coord_equal() +
    scale_fill_manual(values=c("black", "white"), guide="none")+
    scale_size_continuous(range=c(5,24),breaks=count.breaks)


}




