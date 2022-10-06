# Server file for Sedeghat insulin model 

library(deSolve)
library(ggplot2)

# First order kinetics

# ui.R

ui <- fluidPage(
  titlePanel("Glucose - insulin model"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Modeling of glucose and insulin concentrations"),
      radioButtons("selection",
                   label   = "Glucose input",
                   choices = c("Intravenous" = 1,
                               "Meal" = 2),
                   selected = 1
      ),
      
      sliderInput("p1", "P1", min = 0.0, max=0.1, value=0.013),
      sliderInput("p2", "P2", min = 0.0, max=0.02, value=0.00347),
      sliderInput("p3", "P3", min = 0.0, max=0.3, value=0.2),  
      sliderInput("p4", "P4", min = 0.0, max=0.1, value=0.03),
     
      
      actionButton("reset", label="Reset")
    ),
    
    mainPanel(
      plotOutput("plot1", height="400px"),
      plotOutput("plot2", height="400px")
    )
  )
)


server <- function(input, output, session) {
    
    # Reset paramater values
    observeEvent(input$reset, { updateSliderInput(session, "p1", value = 0.013)
                                updateSliderInput(session, "p2", value = 0.00347)
                                updateSliderInput(session, "p3", value = 0.2)
                                updateSliderInput(session, "p4", value = 0.03)
                                })
  
  
  bolie_iv <- function (t, y, p) {
    dose    <- 75
    t_dose  <- 20
    expon   <- (t - t_dose)^2 / 4
    R       <- dose * exp(-expon)
    dgdt    <- -p[3]*y[2] - p[4]*y[1] + R  # glucose
    didt    <- -p[1]*y[2] + p[2]*y[1]    # insulin 
 
     return(list(c(dgdt, didt)))
  }
  
  bolie_meal <- function(t, y, p) {
    # Solves ODE's for meal Bolie model
    a        <- 0.2; 
    b        <- 0.4;
    D        <- 400;       # mg/dl assuming 17.5 l plasma volume
    t_dose   <- 20;   
    
    if (t > t_dose) {
      R = D*((b-a)*(b+a)/(2*a))*exp(-0.5*b*(t - t_dose))*sinh(0.5*a*(t - t_dose))
    }
    else {
      R = 0
    }
    dgdt = -p[3]*y[2] - p[4]*y[1] + R  # glucose
    didt = -p[1]*y[2] + p[2]*y[1]      # insulin
    return(list(c(dgdt, didt)))
  }                           
    
      bolie_model <- reactive ({
        # Calculates insulin and glucose concentrations with the 
        # Bolie model
        
        t  <- seq(0, 200, by = .05)
        y0 <- c(300, 0)
        p  <- c(input$p1, input$p2, input$p3, input$p4)
        select <- input$selection
        
         if(select == 1) {                   # intravenous
           out <- ode(y = y0, times=t, func=bolie_iv, parms = p)
         }
         else if(select == 2) {              # meal
           out <- ode(y = y0, times=t, func=bolie_meal, parms = p)
        }
        
        #Results are difference with basal levels, hence add basal
        # concentration: 80 mg/ml for glucose, 10 mu U/ml for insulin
       
        G_corr <- 80.0 + out[, 2]
        I_corr <- out[, 3]
        glucose.df <- cbind.data.frame(out[, 1], G_corr)
        insulin.df <- cbind.data.frame(out[, 1], I_corr)
        
      
        p <- ggplot(data = glucose.df, aes(x = glucose.df[ ,1],
                                           y = glucose.df[ ,2]) )
        p <- p + geom_line(col = "blue", lwd = 2)
        
        p <- p + labs(title = "Glucose concentration",  
                           x = 'Time (minutes)', 
                           y = 'Concentration (mg/dL)') 
        
        p <- p +   theme(text = element_text(size = 14),
                     legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 16),
                     plot.title = element_text(size = 16, 
                     face="bold", hjust = 0.5)
                   )
        
        
        q <- ggplot(data = insulin.df, aes(x = insulin.df[, 1],
                                           y = insulin.df[, 2]))
        q <- q + geom_line(col = "blue", lwd = 2)
        
        q <- q + labs(title = "Insulin concentration",  
                      x = 'Time (minutes)', y = 'Concentration (micro U/L)') 
        q <- q +   theme(text = element_text(size = 14),
                     legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 16),
                     plot.title = element_text(size = 16, 
                     face="bold", hjust = 0.5)
                  )

 
      return(list(p, q))
    })
    
    output$plot1 <- renderPlot(print(bolie_model() [[1]] ))
    output$plot2 <- renderPlot(print(bolie_model() [[2]] ))
  }


shinyApp(ui, server)