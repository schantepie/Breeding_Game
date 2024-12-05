##  The Breeding Game is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License v.3.0 for more details.



library(shiny)
library(DT)
library(ggplot2)
library(kinship2)
library(gridExtra)
##############################################
####### Own R function for the breeding game
##############################################
maxInbreedCoef=0.3
mean_repro_allowed=12
set.seed(666)


c100 <- c(rep(
  c("dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"),times=4)
)




estimatePhenotype<-function(data,generation,nb_loci,varianceEnvironmental,meantrait){
  Generation=rep(generation,(dim(data)[1]))
  BreedingValue=apply(data[,1:(nb_loci*2),drop=FALSE],1,sum)
  EnvironmentalValue= rnorm(dim(data)[1],mean=0,sd=sqrt(varianceEnvironmental))
  PhenotypicValue=BreedingValue+EnvironmentalValue
  if(!is.null(meantrait)){
    PhenotypicValue=PhenotypicValue+rep(meantrait,dim(data)[1])
  }
  Remaininfoffsprigbeforedead=rpois(dim(data)[1],mean_repro_allowed)
  data=cbind(data,Generation,BreedingValue,EnvironmentalValue,PhenotypicValue,Remaininfoffsprigbeforedead,tempOffproduced=rep(0,dim(data)[1]))
  if(generation==0) rownames(data)=paste("Ind_gen_0",1:(dim(data)[1]),sep="_")
  return(data)
}


# classic QTL reproduction
repro<-function(dataset,reproducer,generationNumber,lastknownOffspring,nb_loci,varEnvironmentale,recombination){
  nb_reproducer=dim(reproducer)[1]
  # output$coupleout<-  renderPrint({nb_reproducer})
  offspring= matrix(NA,nrow=nb_reproducer,ncol=nb_loci*2)
  seqFirstAllele=1:nb_loci
  seqSecondAllele=(nb_loci+1):(nb_loci*2)
  
  if(recombination==0.5){
    
    locisampled=matrix(rbinom(nb_reproducer*nb_loci*2,1,0.5),nrow=nb_reproducer)
    offspring[,seqFirstAllele]=dataset[match(reproducer[,1],rownames(dataset)),seqFirstAllele]*locisampled[,seqFirstAllele]+ 
      dataset[match(reproducer[,1],rownames(dataset)),seqSecondAllele]*(+(!locisampled[,seqFirstAllele]))
    offspring[,seqSecondAllele]=dataset[match(reproducer[,2],rownames(dataset)),seqFirstAllele]*locisampled[,seqSecondAllele]+
      dataset[match(reproducer[,2],rownames(dataset)),seqSecondAllele]*(+(!locisampled[,seqSecondAllele]))
  }else{ ##choose among chromosome
    chromosomesampled=matrix(rbinom(nb_reproducer*2,1,0.5),nrow=nb_reproducer)
    chromosomesampled=matrix(chromosomesampled[rep(1:nb_reproducer, 2), ], nrow =nb_reproducer, ncol =nb_loci*2)
    offspring[,seqFirstAllele]=dataset[match(reproducer[,1],rownames(dataset)),seqFirstAllele]*chromosomesampled[,seqFirstAllele]+
      dataset[match(reproducer[,1],rownames(dataset)),seqSecondAllele]*(+(!chromosomesampled[,seqFirstAllele]))
    offspring[,seqSecondAllele]=dataset[match(reproducer[,2],rownames(dataset)),seqFirstAllele]*chromosomesampled[,seqSecondAllele]+
      dataset[match(reproducer[,2],rownames(dataset)),seqSecondAllele]*(+(!chromosomesampled[,seqSecondAllele]))
  }
  
  offspring=estimatePhenotype(offspring,generation=generationNumber+1,nb_loci = nb_loci,varianceEnvironmental = varEnvironmentale , meantrait=NULL)
  # offspring=cbind(offspring,InbreedCoef=0,alive=rep(1,dim(offspring)[1]),deadAtGeneration=NA)
  # offspring=liveOrDead(offspring,generationNumber)
  rownames(offspring)=paste("Ind_gen",generationNumber+1,lastknownOffspring():(lastknownOffspring()+(dim(offspring)[1])-1),sep="_")
  lastknownOffspring(lastknownOffspring()+(dim(offspring)[1]))
  return(offspring)
}



################################
####### Founders initialization
################################
founderSet_UI <- function(id){
  ns <- NS(id)
  tagList(
    h1("Generate the population of founders"),
    sliderInput(ns("meantrait"), "Choose the mean trait at start", min=1, max=5, value = c(2),step=1),
    sliderInput(ns("nb_loci"), "Choose the number of loci", min=1, max=30, value = c(10),step=1),
    sliderInput(ns("varianceEnvironmentale"), "Choose an environmental variance", min=1, max=5, value = c(2),step=0.1),
    sliderInput(ns("varianceGenetic"), "Choose a genetic variance", min=1, max=5, value = c(2),step=0.1),
    sliderInput(ns("nbFounders"), "Choose a number of founders", min=1, max=100, value = c(50),step=10),
    # sliderInput(ns("nbFounders"), "Choose a number of founders", min=1, max=100, value = c(50),step=10),
    # sliderInput(ns("recombination"), "Choose a recombination rate:", min=0, max=1, value = c(0),step=0.5),
    
    actionButton(ns(paste0("founderpop")), "Initialize/Renitialize the population of founders"),
    uiOutput(ns("datareact")),
    uiOutput(ns("pedi")),
    verbatimTextOutput(ns("ula")),
    h3(htmlOutput(ns("textfounders"))),
    h3(htmlOutput(ns("textfounders2"))),
    actionButton(ns(paste0("addgenerationfounder")), "Start breeding")
  )
}

founderInitialisation_server <- function(id,datareact,current_generation,number_loci,varEnvironmentale,parent,buttongeneration){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    observeEvent(input[[paste0("founderpop")]],{
      ###### Some value used outsite of the module
      number_loci(input$nb_loci)
      varEnvironmentale(input$varianceEnvironmentale)
      ########################################
      if(current_generation()==-1){
        ##### deal with founder population
        data_founder=matrix(rnorm(2*input$nb_loci*input$nbFounders,mean=0, sd=sqrt(input$varianceGenetic/(input$nb_loci))),nrow=input$nbFounders)
        data_founder=estimatePhenotype(data_founder,generation=current_generation()+1,nb_loci = input$nb_loci,varianceEnvironmental = input$varianceEnvironmentale,meantrait=input$meantrait)
        datareact[["allGenerations"]]<-data_founder
        datareact[["pedigree"]]<-cbind(ID=rownames(data_founder),Mother=NA,Father=NA,sex=sample(c(1,2),size = dim(data_founder)[1],replace = T),InbreedCoef=0,alive=rep(1,dim(data_founder)[1]),generation=0)
        
        
        # output$datareact<- renderTable(datareact$allGenerations,rownames=TRUE)
        # output$datareact<- renderTable(datareact$tempOffsprings,rownames=TRUE)
        # output$pedi<- renderTable(datareact[["pedigree"]][,c("generation","InbreedCoef")],rownames=TRUE)
        # output$pedi<- renderTable( datareact[["tidyp"]])

        
        if(input[[paste0("founderpop")]]==1){
          output$textfounders <- renderUI({HTML(paste("The founders are initialized."," ", sep="<br/>")) })
          output$textfounders2 <- renderUI({HTML(paste(" ", sep="<br/>")) })
        }else{
          output$textfounders <- renderUI({HTML(paste(paste0("The founders are reinitialized : #",input$founderpop-1)," ", sep="<br/>")) })
          output$textfounders2 <- renderUI({HTML(paste(" ", sep="<br/>")) })
        }
        
      }else{
        showNotification("NOPE, to reinialize founders, remove all the generations previously created", type = "error",closeButton = TRUE)
        output$textfounders <- renderUI({HTML(paste("The founders are not reinitialized", sep="<br/>")) })
      }
    })
    
    observeEvent(input[[paste0("addgenerationfounder")]], {
      if(is.null(datareact[["allGenerations"]])){
        output$textfounders2 <- renderUI({HTML(paste("Hummm!! Initialize founders first !!!"," ", sep="<br/>")) })
      }else{
        if(current_generation()==-1){
          buttongeneration(rnorm(1))
          
        }else{
          showNotification("NOPE, to restart the game, remove all the generations previously created", type = "error",closeButton = TRUE)
        }
      }
      
    })
    
    
    
  })
}


################################
####### Select pairs to reproduce
################################

couple_input_UI <- function(id,datareact,cpt){
  ns <- NS(id)
  # namesInd=rownames(datareact$allGenerations)[which(datareact$allGenerations[,"alive"]==1)]
  namesMale = datareact$pedigree[which(datareact$pedigree[,"alive"]==1 & datareact$pedigree[,"sex"]==1),"ID"]
  namesFemale = datareact$pedigree[which(datareact$pedigree[,"alive"]==1 & datareact$pedigree[,"sex"]==2),"ID"]
  tagList(
    h1("1. Selecting pairs for reproduction"),
    fluidRow( column(2,selectizeInput(inputId=ns('selectedFemale'),
                                      label='Select a female to form a pairs',
                                      choices = namesFemale,
                                      selected = NULL,
                                      multiple= TRUE , 
                                      options = list(minItems = 1,maxItems = 1))),
              column(2,selectizeInput(inputId=ns('selectedMale'),
                                      label='Select a male to form a pairs',
                                      choices = namesMale,
                                      selected = NULL,
                                      multiple= TRUE ,
                                      options = list(minItems = 1,maxItems = 1)),actionButton(ns(paste0("combine",cpt)), "Add a couple"))
              , 
              
              # column(3,uiOutput(ns("del"))),
              column(5, checkboxGroupInput(ns("select"), "Select pairs to delete" ),
                     actionButton(ns(paste0("dodelete",cpt)), HTML("Delete selected <br/> click without selection to remove all pairs",sep="<br/>")))
                   
)
  )
}

couple_input_server <- function(id,clean,cpt,datareact,couplereact,current_generation,generationFixed,hide=FALSE){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    if(clean==FALSE){
      
      observeEvent(input[[paste0("combine",cpt)]], {
   
        
        if(current_generation()==generationFixed){
          
          reproF=reproM=FALSE
          if(is.null(input$selectedFemale)){
            showNotification("You should select a Female",type="error")
          }else if (datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedFemale),"Remaininfoffsprigbeforedead"]== datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedFemale),"tempOffproduced"]){
            showNotification(paste0(input$selectedFemale, " cannot reproduce anymore. She already produced too much offsprings and died of love"),type="error")
          }else{
            reproF=TRUE
          }
          if(is.null(input$selectedMale)){
            showNotification("You should select a Male",type="error")
          }else if (datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedMale),"Remaininfoffsprigbeforedead"]== datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedMale),"tempOffproduced"]){
            showNotification(paste0(input$selectedMale, " cannot reproduce anymore. He already produced too much offsprings and died of love"),type="error")
          } else{
            reproM=TRUE
          }
          
          if(sum(reproF,reproM)==2){
            datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedFemale),"tempOffproduced"]=
              datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedFemale),"tempOffproduced"]+1
            datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedMale),"tempOffproduced"]=
              datareact[["allGenerations"]][which(rownames(datareact$allGenerations)==input$selectedMale),"tempOffproduced"]+1
            
            couplereact[["couple"]]<-append(couplereact[["couple"]], list(c(input$selectedFemale,input$selectedMale)))
            couplereact[["coupletable"]]<-do.call(rbind,couplereact[["couple"]])
            couplereact[["Names"]]<-paste(seq_along( couplereact[["couple"]]),": (" ,  couplereact[["coupletable"]][,1]," & ",  couplereact[["coupletable"]][,2],")",sep="")
          }  
        }else{
          showNotification("NOPE, remove all the next generations before",type="error")
        }
      })
      
      
      observeEvent(input[[paste0("dodelete",cpt)]],{
        if(!is.null( couplereact[["couple"]])){
          if(is.null(input$select)){
            # datareact[["allGenerations"]][match(unique(c(couplereact[["coupletable"]])),rownames(datareact$allGenerations)),"tempOffproduced"]= 0
            freqremove=as.data.frame(table(couplereact[["coupletable"]]))
            datareact[["allGenerations"]][match(freqremove[,1],rownames(datareact$allGenerations)),"tempOffproduced"]=  datareact[["allGenerations"]][match(freqremove[,1],rownames(datareact$allGenerations)),"tempOffproduced"]-freqremove[,2]
            couplereact[["couple"]]=NULL
            couplereact[["coupletable"]]=NULL
            couplereact[["Names"]]="No pairs remaining"
          }else{
            removeselect=as.numeric(do.call(rbind,strsplit(input$select, split = ":"))[,1])
            freqremove=as.data.frame(table(couplereact[["coupletable"]][as.integer(removeselect),]))
            datareact[["allGenerations"]][match(freqremove[,1],rownames(datareact$allGenerations)),"tempOffproduced"]=  datareact[["allGenerations"]][match(freqremove[,1],rownames(datareact$allGenerations)),"tempOffproduced"]-freqremove[,2]
            couplereact[["couple"]] = couplereact[["couple"]][-as.integer(removeselect)]
            couplereact[["coupletable"]]=do.call(rbind,couplereact[["couple"]] )
            couplereact[["Names"]]=paste(seq_along(couplereact[["couple"]] ),": (" , couplereact[["coupletable"]][,1]," & ", couplereact[["coupletable"]][,2],")",sep="")
            if(length( couplereact[["couple"]])==0) couplereact[["Names"]]="No pairs remaining"
          }
        }
      })
      observe({
        updateCheckboxGroupInput(session, "select",
                                 label = NULL,
                                 choices =  couplereact[["Names"]],
                                 selected = NULL)
      })
      return(couplereact[["coupletable"]])
    }else{
      couplereact[["couple"]]=NULL
      couplereact[["coupletable"]]=NULL
      couplereact[["Names"]]="No pairs remaining"
      datareact[["allGenerations"]][match(unique(c(couplereact[["coupletable"]])),rownames(datareact$allGenerations)),"tempOffproduced"]= 0
    }
    
    
  })
}




################################
####### Reproduction of Pairs
################################

reproduction_UI <- function(id,cpt){
  ns <- NS(id)
  tagList(
    h1("2. Perform reproduction"),
    # sliderInput(ns("recombination"), "Choose a recombination rate:", min=0, max=1, value = c(0),step=0.5),
    actionButton(ns(paste0("repro",cpt)), "Reproduce the pairs"),
    actionButton(ns(paste0("resetrepro",cpt)), "Kill all the produced offsprings"),
    verbatimTextOutput(ns("reproducer")),
    h3(htmlOutput(ns("numberOffinQueue"))),
    uiOutput(ns("offPlot")),
    verbatimTextOutput(ns("coupleout"))
  )
}




reproduction_server <- function(id,current_generation,datareact,nb_loci,varEnvironmentale,newOff_dimensions,cpt,couplereact){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    lastknownOffspring <- reactiveVal(1)

  
    observeEvent(input[[paste0("repro",cpt)]], {
      if(!is.null(couplereact[["coupletable"]])){
        # message(input$recombination)
        newOff<-repro(dataset = datareact$allGenerations, reproducer = couplereact[["coupletable"]],
                      generationNumber = current_generation,lastknownOffspring = lastknownOffspring,
                      nb_loci = nb_loci,varEnvironmentale=varEnvironmentale(),recombination=0.5)
        
        datareact[["tempOffsprings"]] <- rbind( datareact[["tempOffsprings"]] ,newOff)

        datareact[["tempPedigree"]] <- rbind(datareact[["tempPedigree"]] ,cbind(ID=rownames(newOff),
                                                                                Mother=couplereact[["coupletable"]][,1],
                                                                                Father=couplereact[["coupletable"]][,2],
                                                                                sex=sample(c(1,2),size = dim(couplereact[["coupletable"]])[1],replace = T),
                                                                                InbreedCoef= datareact[["longkin"]][match( paste0(couplereact[["coupletable"]][,1],couplereact[["coupletable"]][,2]),
                                                                                                                           paste0(datareact[["longkin"]][,"Mother"],datareact[["longkin"]][,"Father"])),"Inbreeding coefficient of produced offsprings"],
                                                                                alive=rep(1,dim(couplereact[["coupletable"]])[1]),
                                                                                generation = current_generation+1)
                                                 )
  
        
        deadOffspringLine=which(datareact[["tempPedigree"]][,"InbreedCoef"]>=maxInbreedCoef)
        if(sum(deadOffspringLine)>0){
          datareact[["tempPedigree"]][datareact[["tempPedigree"]][,"InbreedCoef"]>=maxInbreedCoef,"alive"]=0
          for(i in datareact[["tempPedigree"]][datareact[["tempPedigree"]][,"alive"]==0,"ID"]){
            showNotification(paste("Offspring",i,sep=" "), " died because he was too inbreed",type = "error",duration = 15)
            datareact[["tempOffsprings"]] <- datareact[["tempOffsprings"]][-which(rownames(datareact[["tempOffsprings"]])==i),,drop=FALSE]
          }
          
          
           datareact[["tempPedigree"]]=datareact[["tempPedigree"]][datareact[["tempPedigree"]][,"InbreedCoef"]<maxInbreedCoef,]
       
          
        }
     
        couplereact[["couple"]]=NULL
        couplereact[["coupletable"]]=NULL
        couplereact[["Names"]]="No pairs remaining"      

 
        
      }
      
    })
    
    
    observeEvent(input[[paste0("resetrepro",cpt)]],{
      datareact[["tempOffsprings"]]=NULL
      datareact[["tempPedigree"]]=NULL
      datareact[["allGenerations"]][,"tempOffproduced"]=0
      lastknownOffspring(1)
      # newOff_dimensions(0)
    })
    
    output$offspringPlot <- renderPlot({        
      bins <- seq(min(datareact[["tempOffsprings"]][,"PhenotypicValue"]),
                  max(datareact[["tempOffsprings"]][,"PhenotypicValue"]), length.out = 10)
      hist(datareact[["tempOffsprings"]][,"PhenotypicValue"], breaks = bins, col = "#007bc2", border = "white",
           xlab = "Phenotypic value",
           main = "Histogram of Phenotypic value for the produced offsprings")
    })
    
    output$needMore <- renderUI({        
      HTML(paste("At least 2 offsprings are needed to display the histogram of phenotypes" , sep="<br/>")) 
    })
    
    output$offPlot <- renderUI({
      if(newOff_dimensions()<2){
        h3(htmlOutput(ns("needMore")))
      }else{
        plotOutput(ns("offspringPlot"))
      }
    })
    output$numberOffinQueue <- renderUI({HTML(paste(paste0(newOff_dimensions()," offprings produced")," ", sep="<br/>")) })
  }
  )
}


################################
####### Space for survival implementation (but to be done)
################################

plotAliveAtGeneration_UI<- function(id){
  ns <- NS(id)
  fluidRow(column(6, plotOutput(ns("distPlot"))),
           column(6, DT::dataTableOutput(ns("PhenoTable")))
  )
}

plotAliveAtGeneration_server<- function(id,datareact,current_generation){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    reducedDatareact    <- data.frame(datareact$allGenerations[which(datareact$allGenerations[,"Generation"]<=current_generation),])
    sexChar=datareact$pedigree[which(datareact$pedigree[,"alive"]==1),"sex"]
    sexChar[sexChar==1]="Male"
    sexChar[sexChar==2]="Female"
    InbreedCoef=datareact$pedigree[which(datareact$pedigree[,"alive"]==1),"InbreedCoef"]
    Family=paste(datareact$pedigree[which(datareact$pedigree[,"alive"]==1),"Mother"],datareact$pedigree[which(datareact$pedigree[,"alive"]==1),"Father"],sep="_")
    Family[Family=="NA_NA"]="Founders"
    reducedDatareact    <- data.frame(reducedDatareact[match(datareact$pedigree[which(datareact$pedigree[,"alive"]==1),"ID"], rownames(reducedDatareact)),
                                                  c("PhenotypicValue"),drop=F],sex=sexChar, family=Family,InbreedCoef=InbreedCoef)


    output$distPlot <- renderPlot({
      bins <- seq(min(reducedDatareact[,"PhenotypicValue"]), max(reducedDatareact[,"PhenotypicValue"]), length.out = 10)
      hist(reducedDatareact[,"PhenotypicValue"], breaks = bins, col = "#007bc2", border = "white",
           xlab = "Phenotypic value",
           main = "Histogram of Phenotypic value")
    })
    output$PhenoTable = DT::renderDataTable({
      datatable(reducedDatareact,
                options = list(
                order = list(list(3, 'family'),list(1, 'PhenotypicValue')))) %>% formatRound(c(1), 2)
    })
  })
}




############################
####### add new generation
############################

addgeneration_UI<-function(id,cpt){
  ns <- NS(id)
  tagList( actionButton(ns(paste0("addgeneration",cpt)), "Go to next generation"),
           h3(htmlOutput(ns("removepriorgen")))
  )
}

addgeneration_server <- function(id,datareact,NewIndproduced,buttongeneration,current_generation,generationFixed,newOff_dimensions,cpt,couplereact){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    observeEvent(input[[paste0("addgeneration",cpt)]], {
      if(current_generation()>as.numeric(generationFixed)){
        showNotification("NOPE, remove all the next generations before",type="error")
      }else if (newOff_dimensions()==0){
        showNotification("At least one offspring is needed to go the next generation",type = "error")
      }else if (sum(as.numeric(datareact[["tempPedigree"]][,"alive"]))==0 ){
        showNotification("All the produced offsprings died because of indreeding, try again",type = "error")
      }else{

        datareact[["allGenerations"]] <- rbind(datareact[["allGenerations"]] ,datareact[["tempOffsprings"]] )
        datareact[["pedigree"]] <- rbind(datareact[["pedigree"]] ,datareact[["tempPedigree"]] )
        
        tableremove=as.data.frame(table(c(datareact[["tempPedigree"]][,c("Mother","Father")])))
        
        datareact[["allGenerations"]][match(tableremove[,1],rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]=
          datareact[["allGenerations"]][match(tableremove[,1],rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]-
          tableremove[,2]
        

        #clean before jumping to a new generation
        couplereact[["couple"]]=NULL
        couplereact[["coupletable"]]=NULL
        couplereact[["Names"]]="No pairs selected"   
        datareact[["tempOffsprings"]] = NULL
        datareact[["tempPedigree"]]= NULL
        datareact[["allGenerations"]][,"tempOffproduced"]=0
        
        buttongeneration(rnorm(1))
        
      }
    })
  })
}

############################
####### update offsprings after generation removal
############################

update_offspring_count_server<- function(id,current_generation,datareact,decreaseCount){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    ind_produced=rownames(datareact[["allGenerations"]][which(datareact[["allGenerations"]][,"Generation"]==current_generation),,drop=FALSE])
    
    tablerecover=as.data.frame(table(c(datareact[["pedigree"]][match(ind_produced,datareact[["pedigree"]][,"ID"]),c("Mother","Father")])))
    # 
    
    message(ind_produced)

    
    if(dim(tablerecover)[1]==0){
      
      datareact[["allGenerations"]][match(ind_produced,rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]=
        datareact[["allGenerations"]][match(ind_produced,rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]+
        1
      
      datareact[["pedigree"]]=datareact[["pedigree"]][-match(ind_produced,datareact[["pedigree"]][,"ID"]),]
      
    }else{
      
      datareact[["allGenerations"]][match(tablerecover[,1],rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]=
      datareact[["allGenerations"]][match(tablerecover[,1],rownames(datareact[["allGenerations"]])),"Remaininfoffsprigbeforedead"]+
      tablerecover[,2]
      
      datareact[["pedigree"]]=datareact[["pedigree"]][-match(ind_produced,datareact[["pedigree"]][,"ID"]),]
    }
  

  })
}


############### pedigree_panel


pedigree_UI <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("pedPlot")),
    DT::dataTableOutput(ns("pedTable"))
    # tableOutput(ns("tableSummaryOutput"))
  )
}

 pedigree_server<- function(id,datareact){
      moduleServer(id, function(input, output, session){

        ns <- session$ns
        
        observeEvent(datareact[["pedigree"]],{

          datareact[["pedigree_Kinship2"]]=pedigree(id = datareact[["pedigree"]][,"ID"], dadid = datareact[["pedigree"]][,"Father"], momid = datareact[["pedigree"]][,"Mother"], sex = as.numeric(datareact[["pedigree"]][,"sex"]))
          datareact[["kinship"]]=kinship(datareact[["pedigree_Kinship2"]])
          longkin=cbind(expand.grid(rownames( datareact[["kinship"]]),rownames( datareact[["kinship"]])),c( datareact[["kinship"]]),c(lower.tri( datareact[["kinship"]])))
          colnames(longkin)=c("parent1","parent2","Inbreeding coefficient of produced offsprings","a")     
          longkin=longkin[longkin[,"a"],]
          sex1=datareact[["pedigree"]][match(longkin[,"parent1"],datareact[["pedigree"]][,"ID"]),"sex"]
          sex2=datareact[["pedigree"]][match(longkin[,"parent2"],datareact[["pedigree"]][,"ID"]),"sex"]
          longkin=cbind(longkin,sex1,sex2)
          longkin=longkin[sex1!=sex2,]
          longkin$Father=longkin$parent1
          longkin[longkin$sex2==1,"Father"]=longkin[longkin$sex2==1,"parent2"]
          longkin$Mother=longkin$parent2
          longkin[longkin$sex1==2,"Mother"]=longkin[longkin$sex1==2,"parent1"]
          
          datareact[["longkin"]]=longkin[,c("Mother","Father","Inbreeding coefficient of produced offsprings")]
          
        #   if(max(datareact[["pedigree"]][,"generation"])>0){
        #       datareact[["tidyp"]]=tidyped(datareact[["pedigree"]][,c("ID","Father","Mother","generation")],addgen=TRUE)
        #       datareact$tidyp$Gen=as.numeric(datareact$tidyp$generation)+1
        # 
        #       
        # }
          
          output$pedPlotvispred <- renderPlot({
            
            plot.pedigree(datareact[["pedigree_Kinship2"]],
                          id=datareact[["pedigree"]][,c("generation")],
                          col=c100[as.numeric(datareact[["pedigree"]][,c("generation")])+1],
                          affected=rep(1,length(datareact[["pedigree"]][,c("generation")])))
           # visped( datareact[["tidyp"]],compact = TRUE)
          })
          
          # message(c100[as.numeric(as.character(datareact[["pedigree"]][,c("generation")]))+1])
          
          output$needMoreGeneration <- renderUI({        
            HTML(paste("At least 1 generation is needed to display the pedigree" , sep="<br/>")) 
          })
          
          output$pedPlot <- renderUI({
            if(max(datareact[["allGenerations"]][,"Generation"])==0){
              h3(htmlOutput(ns("needMoreGeneration")))
            }else{
              plotOutput(ns("pedPlotvispred"))
            }
          })
          
          
          
          
          
          
             output$pedTable = DT::renderDataTable({
               datatable(datareact[["longkin"]])
             })


        })
      })
 }

########################
####### Results
#######################

Results_UI <- function(id){
  ns <- NS(id)
  tagList(
    plotOutput(ns("violonplotpheno")),
    plotOutput(ns("violonplotVa")),
    plotOutput(ns("violonplotInbreeding")),
    tableOutput(ns("tableSummaryOutput")),
    tableOutput(ns("gener")),
    downloadButton(ns('downloadResults'), 'Download the results')
  )
}

results_server <- function(id,datareact){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    summaryStat<-reactiveValues()
    exportresults <- reactiveValues(plotpheno=NULL,
                            plotgenet=NULL,
                            plotinbreed=NULL,
                            tableResults=NULL)
    
    observeEvent(datareact[["allGenerations"]],{
      # (max(as.numeric(Generation))+1)
      output$violonplotpheno<-renderPlot({
        exportresults$plotpheno<-ggplot(as.data.frame(datareact[["allGenerations"]]),aes(x=factor(Generation), y=PhenotypicValue)) + 
          ylab("Phenotypic values")+
          xlab("Generation")+
          geom_violin(aes(fill=factor(Generation)), draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1)+
          scale_fill_manual(values=c100[unique(datareact[["allGenerations"]][,"Generation"])+1])+
          # geom_boxplot(width=0.3)+
          stat_summary(fun=mean, geom="point", shape=20, size=5)+
          ggtitle("Violon plot of phenotypic values across generations") +
          theme(text = element_text(size = 20),legend.position="none")
        exportresults$plotpheno
      })
      
    
      
      output$violonplotVa<-renderPlot({
        exportresults$plotgenet<-ggplot(as.data.frame(datareact[["allGenerations"]]),aes(x=factor(Generation), y=BreedingValue)) + 
          ylab("Breeding values")+
          xlab("Generation")+
          geom_violin(aes(fill=factor(Generation)), draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1)+
          scale_fill_manual(values=c100[unique(datareact[["allGenerations"]][,"Generation"])+1])+
          # geom_boxplot(width=0.3)+
          stat_summary(fun=mean, geom="point", shape=20, size=5)+
          ggtitle("Violon plot of breeding values across generations") +
          theme(text = element_text(size = 20),legend.position="none")
        exportresults$plotgenet
      })
  
     summaryStat[["table"]]<- cbind(with(data.frame(datareact[["allGenerations"]][,c("Generation","PhenotypicValue","BreedingValue")]),
                                        aggregate(cbind(PhenotypicValue,BreedingValue) , by=list(as.factor(Generation)) , mean)),
                                     with(data.frame(datareact[["allGenerations"]][,c("Generation","PhenotypicValue","BreedingValue")]),
                                        aggregate(cbind(PhenotypicValue,BreedingValue) , by=list(as.factor(Generation)) , var)) ,
                                     with(data.frame(datareact[["allGenerations"]][,c("Generation","PhenotypicValue")]),
                                          aggregate(cbind(PhenotypicValue) , by=list(as.factor(Generation)) , function(x){length(x)}))
                             )
       colnames(summaryStat$table)<-c("Generation","Mean Phenotype","Mean Breeding Value","Generation2","Phenotypic variance", "Genetic variance","Generation3","Nb individuals")
       summaryStat$table=summaryStat$table[,c("Generation","Mean Phenotype","Phenotypic variance","Mean Breeding Value","Genetic variance","Nb individuals")]
       exportresults$tabResults<-tableGrob(summaryStat[["table"]])
       
      output$tableSummaryOutput <- renderTable({
        summaryStat[["table"]]
      })
    })
    
    
    observeEvent(datareact[["pedigree"]],{
     
      output$violonplotInbreeding<-renderPlot({
        datareact[["pedigree"]][,"generation"]=as.numeric(datareact[["pedigree"]][,"generation"])
        exportresults$plotinbreed<-ggplot(as.data.frame(datareact[["pedigree"]]),aes(x= reorder(generation, sort(as.numeric(generation))), y=InbreedCoef)) + 
          ylab("Inbreeding coefficients")+
          xlab("Generation")+
          geom_violin(aes(fill=factor(generation)), draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(height = 0, width = 0.1)+
          scale_fill_manual(values=c100[unique(datareact[["allGenerations"]][,"Generation"])+1])+
          # geom_boxplot(width=0.3)+
          stat_summary(fun=mean, geom="point", shape=20, size=5)+
          ggtitle("Violon plot of inbreeding coefficient across generations") +
          # scale_x_continuous(breaks = unique(datareact[["pedigree"]][,"generation"]))+
          theme(text = element_text(size = 20),legend.position="none")
        exportresults$plotinbreed
      })
      
      
    
      
      ## clicking on the export button will generate a pdf file 
      ## containing all stored plots and tables
      output$downloadResults = downloadHandler(
        filename = function() {paste("Results-", Sys.Date(), ".pdf", sep="")},
        content = function(file) {
          
          pdf(file, onefile = TRUE, width = 19, height = 25)
          grid.arrange(exportresults$plotpheno,
                       exportresults$plotgenet,
                       exportresults$plotinbreed,
                       exportresults$tabResults,
                       nrow = 4,
                       ncol = 1)
          plot.pedigree(datareact[["pedigree_Kinship2"]],
                        id=datareact[["pedigree"]][,c("generation")],
                        col=c100[as.numeric(datareact[["pedigree"]][,c("generation")])+1],
                        affected=rep(1,length(datareact[["pedigree"]][,c("generation")])))
          dev.off()
          
          
        
      # if(max(datareact[["allGenerations"]][,"Generation"])>0){
      #   output$downloadResults = downloadHandler(
      #     filename = function() {paste("Pedigree-", Sys.Date(), ".pdf", sep="")},
      #     content = function(file) {
      #        visped( datareact[["tidyp"]],compact = TRUE,file=file)
      #     })
      # }
    })
    
    
    
  })
})
}


startPage_UI <- function(id){
  ns <- NS(id)
  tagList(  h3(htmlOutput("instructions")),
            actionButton("startBreedingGame", "Start the Breeding Game", icon("paper-plane"), 
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}

TabPanelSetUp_UI<- function(id){
  ns <- NS(id)
  actionButton("removelastGen", "Remove the last generation", icon = icon("minus-circle"))
  tabsetPanel(id = "allgenerations",
              tabPanel("Founders",founderSet_UI("founders")),
              tabPanel("Results",Results_UI("results")))
}

ui <- fluidPage(
  startPage_UI("BreedingGame"),
  tabsetPanel(id = "allgenerations")
)


server <- function(input, output, session) {
                  current_generation <- reactiveVal(-1)
                  datareact<-reactiveValues()
                  pedigree<-reactiveValues()
                  couplereact<- reactiveValues()
                  
                  
                  
                  
                  number_loci<-reactiveVal(NULL)
                  varEnvironmentale<-reactiveVal(NULL)
                  buttongeneration <- reactiveVal(NULL)
                  cpt <- reactiveVal(0)
                  
                  #### Display instructions
                  
                  output$instructions<-renderUI({HTML(paste("Objective: In 1h, try to get the max values for the mean phenotype of your last generation",
                                                            "1. The last generation should be represented by a minimum of 50 individuals",
                                                            "2. You can delete generations and retry if your are not satisfied by your trial",
                                                            "3. Once you are finished with a breeding trial, you can save a pdf of your results located at the bottom of the results panel", sep="<br/> <br/>")) })
                  observeEvent(input$startBreedingGame,{
                    removeUI(paste0("#instructions"))
                    removeUI(paste0("#startBreedingGame"))
                    insertTab(inputId = "allgenerations", tabPanel("Founders",founderSet_UI("founders")),
                              select=TRUE,target = NULL)
                    
                  })
                  
                  
                  #### Initialize founder population
                  founderInitialisation_server("founders",datareact,current_generation,number_loci,varEnvironmentale,parent=session,buttongeneration)
                  results_server("results",datareact)
                  pedigree_server("pedigree",datareact)
                   
                  observeEvent(buttongeneration(), {
                    insertUI("#allgenerations",ui=actionButton("removelastGen", "Remove the last generation", icon = icon("minus-circle")),  where = c("beforeBegin"))
                    insertTab(inputId = "allgenerations",tabPanel("Results",Results_UI("results")),
                              select=TRUE,target = NULL)
                    insertTab(inputId = "allgenerations",tabPanel("Pedigree",pedigree_UI("pedigree")),
                              select=FALSE,target = "Results",position = "before")
                   
                  },once=TRUE)
                  
                  #### generate the generations 
                  
                  observeEvent(buttongeneration(), {
                    
                    current_generation(current_generation()+1)
                    generationFixed <- sprintf('%s',current_generation())
                    id <- sprintf('generation_%s', current_generation())
                    
                    newOff_dimensions  <- reactive(sum(datareact[["tempOffsprings"]][,"Generation"]==(current_generation()+1)))
                    
                    cpt(cpt()+1)
                    plotAliveAtGeneration_server(id,datareact,current_generation())
                    
                    insertTab(inputId = "allgenerations",
                              tab=tabPanel(id,
                                           plotAliveAtGeneration_UI(id),
                                           couple_input_UI(id, datareact,cpt()),
                                           reproduction_UI(id,cpt()),
                                           addgeneration_UI(id,cpt())
                              ),
                              select=TRUE,
                              target="Results",
                              position = c("before")
                    )
                    # 
                    couple_input_server(id,clean=FALSE,cpt(),datareact,couplereact,current_generation,generationFixed)
                    NewIndproduced<-reproduction_server(id,current_generation(),datareact,nb_loci = number_loci(),varEnvironmentale=varEnvironmentale,newOff_dimensions,cpt(),couplereact)
                    addgeneration_server(id,datareact,NewIndproduced,buttongeneration,current_generation,generationFixed,newOff_dimensions,cpt(),couplereact)
                    # pedigree_update_server("Pedigree",datareact)
                  })
                  
                  
                  observeEvent(input$removelastGen, {
                    id <- sprintf('generation_%s', current_generation())
                    
                    if (current_generation()>0) {
                      # if (current_generation()>=1) {update_offspring_count_server(id,current_generation(),datareact,decreaseCount = FALSE)}
                      update_offspring_count_server(id,current_generation(),datareact,decreaseCount = FALSE)
                      datareact[["allGenerations"]]<- datareact[["allGenerations"]][which(datareact[["allGenerations"]][,"Generation"]<current_generation()),]
                      couplereact[["couple"]]=NULL
                      couplereact[["coupletable"]]=NULL
                      couplereact[["Names"]]="No pairs selected"   
                      
                    }
                    removeTab(inputId = "allgenerations", target = id )
                    
                    
                    
                    current_generation(current_generation() - 1 )
                    if (current_generation()<(-1)) current_generation(-1)
                  })
                  
                  
                  
                  
  
}

# Run the application 
shinyApp(ui = ui, server = server)




