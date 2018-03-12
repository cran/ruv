ruv_shiny <-
function (Y, rowinfo, colinfo, options = list(port = 3840)) 
{
    if (!(requireNamespace("shiny", quietly = TRUE) & requireNamespace("colourpicker", 
        quietly = TRUE))) 
        stop("The ruv_shiny function requires the following packages: shiny, colourpicker, scales")
    Y.input = Y
    rowinfo.input = rowinfo
    colinfo.input = colinfo
    shiny::shinyApp(ui = shiny::fluidPage(shiny::navbarPage("RUV", 
        shiny::tabPanel("Data", shiny::fluidRow(shiny::column(3, 
            shiny::wellPanel(shiny::fileInput("Y.file", label = "Data Matrix", 
                accept = c("text/csv", "text/comma-separated-values,text/plain", 
                  ".csv")), shiny::fileInput("rowinfo.file", 
                label = "Row Information", accept = c("text/csv", 
                  "text/comma-separated-values,text/plain", ".csv")), 
                shiny::fileInput("colinfo.file", label = "Column Information", 
                  accept = c("text/csv", "text/comma-separated-values,text/plain", 
                    ".csv")))), shiny::column(9, shiny::h2("Data Matrix"), 
            shiny::dataTableOutput("datamatrixtable"), shiny::h2("Row Information"), 
            shiny::dataTableOutput("rowtable"), shiny::h2("Column Information"), 
            shiny::dataTableOutput("coltable")))), shiny::tabPanel("Canonical Correlation", 
            shiny::fluidRow(shiny::column(3, shiny::wellPanel(shiny::uiOutput("Xvar.CC"), 
                shiny::uiOutput("ctl.CC"))), shiny::column(9, 
                shiny::h2("Canonical Correlation Plot"), shiny::plotOutput("ccplot", 
                  height = "1000px")))), shiny::tabPanel("Global Adjustment", 
            shiny::fluidRow(shiny::column(3, shiny::wellPanel(shiny::h2("Analysis Parameters"), 
                shiny::h3("Centering"), shiny::checkboxInput("centercol", 
                  "Center Columns", value = TRUE), shiny::checkboxInput("centerrow", 
                  "Center Rows by Neg. Controls", value = FALSE), 
                shiny::uiOutput("ctl.I"), shiny::h3("RUV-III"), 
                shiny::uiOutput("ctl.III"), shiny::uiOutput("M"), 
                shiny::uiOutput("burst"), shiny::radioButtons("k.III.select", 
                  "Selection of k", choices = c("Maximum Possible", 
                    "Select Manually"), selected = "Select Manually"), 
                shiny::uiOutput("k.III.sel"), shiny::tags$b("Average Replicates?"), 
                shiny::checkboxInput("avg.III", "Average", value = FALSE), 
                shiny::h2("Display Parameters"), shiny::checkboxGroupInput("plot.k", 
                  "SVD Plot k", choices = 1:5, selected = 1:3), 
                shiny::uiOutput("shape.row"), shiny::uiOutput("color.by.row"), 
                shiny::checkboxInput("ccolor.row", "Custom Coloring", 
                  value = FALSE), shiny::uiOutput("customcolor.row"), 
                shiny::uiOutput("color.by.col"), shiny::checkboxInput("ccolor.col", 
                  "Custom Coloring", value = FALSE), shiny::uiOutput("customcolor.col"))), 
                shiny::column(9, shiny::h2("RLE Plots"), shiny::tabsetPanel(shiny::tabPanel("Centered", 
                  shiny::plotOutput("rle1")), shiny::tabPanel("RUVIII", 
                  shiny::plotOutput("rle3"))), shiny::h2("SVD Plot"), 
                  shiny::tabsetPanel(shiny::tabPanel("Centered", 
                    shiny::plotOutput("svd.plot1", height = "1000px")), 
                    shiny::tabPanel("RUVIII (embedded)", shiny::plotOutput("svd.plot31", 
                      height = "1000px")), shiny::tabPanel("RUVIII", 
                      shiny::plotOutput("svd.plot3", height = "1000px")))))), 
        shiny::tabPanel("Differential Expression", shiny::fluidRow(shiny::column(3, 
            shiny::wellPanel(shiny::h2("Analysis Parameters"), 
                shiny::uiOutput("DE.Y.var"), shiny::uiOutput("Xvar"), 
                shiny::uiOutput("ctl"), shiny::radioButtons("DE.method", 
                  "Method", choices = c("Unadjusted", "RUV2", 
                    "RUV4", "RUVinv", "RUVrinv"), selected = "RUVrinv"), 
                shiny::uiOutput("k"), shiny::tags$b("Variance estimates"), 
                shiny::radioButtons("vtype", label = NULL, choices = c("Standard", 
                  "E. Bayes (Limma)", "Pooled"), selected = "E. Bayes (Limma)"), 
                shiny::tags$b("t statistic"), shiny::radioButtons("ptype", 
                  label = NULL, choices = c("Standard", "Rescaled", 
                    "Empirical"), selected = "Standard"), shiny::h2("Display Parameters"), 
                shiny::uiOutput("pvar"), shiny::uiOutput("bvar"), 
                shiny::uiOutput("color.by"), shiny::checkboxInput("ccolor", 
                  "Custom Coloring", value = FALSE), shiny::uiOutput("customcolor"), 
                shiny::uiOutput("infocols"))), shiny::column(9, 
            shiny::h2("P-Value Histogram and ECDF"), shiny::fluidRow(shiny::column(5, 
                shiny::plotOutput("hist")), shiny::column(7, 
                shiny::plotOutput("ecdf"))), shiny::h2("Top-Ranked Positive Controls"), 
            shiny::fluidRow(shiny::column(10, shiny::plotOutput("rank")), 
                shiny::column(2, shiny::wellPanel(shiny::uiOutput("pctl"), 
                  shiny::uiOutput("rankmaxx"), shiny::uiOutput("rankmaxy")))), 
            shiny::h2("Projection Plot"), shiny::fluidRow(shiny::column(10, 
                shiny::plotOutput("pplot", brush = "pplot_brush")), 
                shiny::column(2, shiny::wellPanel(shiny::uiOutput("ppfactor"), 
                  shiny::uiOutput("ppadjust")))), shiny::tableOutput("ppbrush_info"), 
            shiny::h2("Volcano and Variance Plots"), shiny::fluidRow(shiny::column(6, 
                shiny::plotOutput("volcano", brush = "volcano_brush"), 
                shiny::tableOutput("volbrush_info")), shiny::column(6, 
                shiny::plotOutput("varplot", brush = "varplot_brush"), 
                shiny::tableOutput("varbrush_info"))))), shiny::fluidRow(shiny::tags$hr(), 
            shiny::h2("Top Ranked"), shiny::column(2, shiny::wellPanel(shiny::uiOutput("statcols"))), 
            shiny::column(10, shiny::dataTableOutput("table")))))), 
        server = function(input, output) {
            options(shiny.maxRequestSize = 100 * 1024^2)
            output$DE.Y.var = shiny::renderUI({
                possible = "Original"
                if (!is.null(Y1())) 
                  possible = c(possible, "Centered")
                if (!is.null(Y3())) 
                  possible = c(possible, "RUVIII")
                shiny::selectInput("DE.Y.var", "Data Matrix", 
                  choices = possible)
            })
            output$Xvar = shiny::renderUI({
                shiny::selectInput("X", "Factor of Interest (X)", 
                  choices = names(rowinfo()))
            })
            output$Xvar.CC = shiny::renderUI({
                shiny::selectInput("X.CC", "Factor", choices = names(rowinfo()))
            })
            output$M = shiny::renderUI({
                shiny::checkboxGroupInput("M", "Replicates", 
                  choices = names(rowinfo()), selected = FALSE)
            })
            output$burst = shiny::renderUI({
                if (length(input$M) == 0) 
                  return(NULL)
                shiny::checkboxGroupInput("burst", "Burst", choices = colnames(replicate.matrix(rowinfo()[, 
                  input$M])), selected = FALSE)
            })
            output$ctl = shiny::renderUI({
                shiny::selectInput("ctl", "Negative Controls", 
                  choices = c(names(colinfo())[unlist(lapply(colinfo(), 
                    is.logical))], "All"))
            })
            output$ctl.CC = shiny::renderUI({
                shiny::selectInput("ctl.CC", "Set", choices = c(names(colinfo())[unlist(lapply(colinfo(), 
                  is.logical))], "All"))
            })
            output$ctl.I = shiny::renderUI({
                shiny::selectInput("ctl.I", "Negative Controls", 
                  choices = c(names(colinfo())[unlist(lapply(colinfo(), 
                    is.logical))], "All"))
            })
            output$ctl.III = shiny::renderUI({
                shiny::selectInput("ctl.III", "Negative Controls", 
                  choices = c(names(colinfo())[unlist(lapply(colinfo(), 
                    is.logical))], "All"))
            })
            output$k = shiny::renderUI({
                if (input$DE.method == "RUV2" | input$DE.method == 
                  "RUV4") 
                  shiny::numericInput("k", "k", value = 1)
                else return(NULL)
            })
            output$k.III.sel = shiny::renderUI({
                if (input$k.III.select == "Select Manually") 
                  shiny::numericInput("k.III", "k", value = 1)
                else return(NULL)
            })
            output$pvar = shiny::renderUI({
                shiny::selectInput("pvar", "P-Value\n(for Histogram, ECDF, Rank plot)", 
                  choices = (c("F", colnames(quick_fit()$misc$X))))
            })
            output$bvar = shiny::renderUI({
                shiny::selectInput("bvar", "X Column\n(for projection, volcano, variance)", 
                  choices = (colnames(quick_fit()$misc$X)))
            })
            output$color.by = shiny::renderUI({
                shiny::selectInput("color.by", "Color by:", choices = c("Negative Controls", 
                  "None", names(colinfo())))
            })
            output$color.by.row = shiny::renderUI({
                shiny::selectInput("color.by.row", "Color by:", 
                  choices = c("None", names(rowinfo())))
            })
            output$shape.row = shiny::renderUI({
                shiny::selectInput("shape.row", "Plot symbols:", 
                  choices = c("None", names(rowinfo())))
            })
            output$color.by.col = shiny::renderUI({
                shiny::selectInput("color.by.col", "Color by:", 
                  choices = c("None", names(colinfo())))
            })
            output$customcolor = shiny::renderUI({
                if (is.null(input$ccolor)) 
                  return(NULL)
                if (!input$ccolor) 
                  return(NULL)
                colorvar = input$color.by
                if (is.null(colorvar)) 
                  return(NULL)
                if (colorvar == "None") 
                  return(NULL)
                if (colorvar == "Negative Controls") 
                  clevels = c("FALSE", "TRUE")
                else clevels = levels(as.factor(colinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(NULL)
                rval = colourpicker::colourInput(paste0(make.names(clevels[1]), 
                  "_color"), clevels[1], value = alpha((scales::hue_pal())(nc)[1], 
                  0.4), allowTransparent = TRUE)
                if (length(clevels) > 1) 
                  for (i in 2:length(clevels)) rval = shiny::tagList(rval, 
                    colourpicker::colourInput(paste0(make.names(clevels[i]), 
                      "_color"), clevels[i], value = alpha((scales::hue_pal())(nc)[i], 
                      0.7), allowTransparent = TRUE))
                return(rval)
            })
            output$customcolor.row = shiny::renderUI({
                if (is.null(input$ccolor.row)) 
                  return(NULL)
                if (!input$ccolor.row) 
                  return(NULL)
                colorvar = input$color.by.row
                if (is.null(colorvar)) 
                  return(NULL)
                if (colorvar == "None") 
                  return(NULL)
                else clevels = levels(as.factor(rowinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(NULL)
                rval = colourpicker::colourInput(paste0(make.names(clevels[1]), 
                  "_color.row"), clevels[1], value = alpha((scales::hue_pal())(nc)[1], 
                  0.4), allowTransparent = TRUE)
                if (length(clevels) > 1) 
                  for (i in 2:length(clevels)) rval = shiny::tagList(rval, 
                    colourpicker::colourInput(paste0(make.names(clevels[i]), 
                      "_color.row"), clevels[i], value = alpha((scales::hue_pal())(nc)[i], 
                      0.7), allowTransparent = TRUE))
                return(rval)
            })
            output$customcolor.col = shiny::renderUI({
                if (is.null(input$ccolor.col)) 
                  return(NULL)
                if (!input$ccolor.col) 
                  return(NULL)
                colorvar = input$color.by.col
                if (is.null(colorvar)) 
                  return(NULL)
                if (colorvar == "None") 
                  return(NULL)
                if (colorvar == "Negative Controls") 
                  clevels = c("FALSE", "TRUE")
                else clevels = levels(as.factor(colinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(NULL)
                rval = colourpicker::colourInput(paste0(make.names(clevels[1]), 
                  "_color.col"), clevels[1], value = alpha((scales::hue_pal())(nc)[1], 
                  0.4), allowTransparent = TRUE)
                if (length(clevels) > 1) 
                  for (i in 2:length(clevels)) rval = shiny::tagList(rval, 
                    colourpicker::colourInput(paste0(make.names(clevels[i]), 
                      "_color.col"), clevels[i], value = alpha((scales::hue_pal())(nc)[i], 
                      0.7), allowTransparent = TRUE))
                return(rval)
            })
            output$pctl = shiny::renderUI({
                shiny::selectInput("pctl", "Positive Controls", 
                  choices = c("None", names(colinfo())[unlist(lapply(colinfo(), 
                    is.logical))]))
            })
            output$rankmaxx = shiny::renderUI({
                shiny::textInput("rankmaxx", "Max x")
            })
            output$rankmaxy = shiny::renderUI({
                shiny::textInput("rankmaxy", "Max y")
            })
            output$ppfactor = shiny::renderUI({
                if (is.null(fit_summary()$misc$ppalpha)) 
                  return(NULL)
                rval = try(shiny::selectInput("ppfactor", "Factor", 
                  choices = c("gradient", 1:nrow(fit_summary()$misc$ppalpha))))
                if (class(rval) == "try-error") 
                  return(NULL)
                return(rval)
            })
            output$ppadjust = shiny::renderUI({
                if (is.null(input$ppfactor)) 
                  return(NULL)
                if (input$ppfactor != "gradient") 
                  shiny::checkboxInput("ppadjusted", "Adjusted", 
                    value = TRUE)
            })
            output$infocols = shiny::renderUI({
                shiny::checkboxGroupInput("infocols", "Variables to include in tables", 
                  choices = c("rownames", names(colinfo())), 
                  selected = c("rownames", names(colinfo())))
            })
            output$statcols = shiny::renderUI({
                if (is.null(quick_fit())) 
                  return(NULL)
                shiny::checkboxGroupInput("statcols", "Results to include in table", 
                  choices = c(names(quick_fit()$C[, 1:(ncol(quick_fit()$C) - 
                    ncol(colinfo()))])), selected = c("F.p"))
            })
            output$ccplot = shiny::renderPlot(try({
                X = rowinfo()[, which(names(rowinfo()) == input$X.CC)]
                if (input$ctl.CC == "All") 
                  ctl = rep(TRUE, ncol(Y()))
                else ctl = colinfo()[, which(names(colinfo()) == 
                  input$ctl.CC)]
                return(ruv_cancorplot(Y(), X, ctl) + theme(text = element_text(size = 20)))
            }))
            output$rle1 = shiny::renderPlot(try({
                ruv_rle(Y1(), rowinfo = rowinfo()) + recolor.rle()
            }))
            output$rle3 = shiny::renderPlot(try({
                if (input$avg.III & !(input$color.by.row %in% 
                  colnames(rowinfo3()))) 
                  return(ruv_rle(Y3()))
                return(ruv_rle(Y3(), rowinfo = rowinfo3()) + 
                  recolor.rle())
            }))
            output$svd.plot1 = shiny::renderPlot(try({
                if (length(plot.k()) == 2) 
                  ruv_svdplot(Y1(), )
                ruv_svdgridplot(Y1(), Y.space = Y1svd(), rowinfo = rowinfo(), 
                  colinfo = colinfo(), Z = NULL, k = plot.k(), 
                  left.additions = svd.left(), right.additions = recolor.col())
            }))
            output$svd.plot3 = shiny::renderPlot(try({
                ruv_svdgridplot(Y3(), Y.space = Y3svd(), rowinfo = rowinfo3(), 
                  colinfo = colinfo(), Z = NULL, k = plot.k(), 
                  left.additions = svd.left3(), right.additions = recolor.col())
            }))
            output$svd.plot31 = shiny::renderPlot(try({
                if (input$avg.III) 
                  return(NULL)
                ruv_svdgridplot(Y3(), Y.space = Y1svd(), rowinfo = rowinfo3(), 
                  colinfo = colinfo(), Z = NULL, k = plot.k(), 
                  left.additions = svd.left(), right.additions = recolor.col())
            }))
            output$hist = shiny::renderPlot(try({
                if (is.null(input$pvar)) 
                  return(NULL)
                if (input$pvar == "") 
                  return(NULL)
                pvar = input$pvar
                if (pvar == "F") 
                  pvar = "all"
                return(ruv_hist(fit_summary(), X.col = pvar))
            }))
            output$ecdf = shiny::renderPlot(try({
                if (is.null(input$pvar)) 
                  return(NULL)
                if (input$pvar == "") 
                  return(NULL)
                pvar = input$pvar
                if (pvar == "F") 
                  pvar = "all"
                return(ruv_ecdf(fit_summary(), X.col = pvar) + 
                  recolor())
            }))
            output$rank = shiny::renderPlot(try({
                if (is.null(input$pvar)) 
                  return(NULL)
                if (input$pctl == "None" | input$pvar == "") 
                  return(NULL)
                fit = fit_summary()
                pctl = fit$C[, input$pctl]
                pvar = input$pvar
                if (pvar == "F") 
                  pvar = "all"
                if (!is.na(as.numeric(input$rankmaxx))) 
                  maxx = as.numeric(input$rankmaxx)
                else maxx = ncol(Y())
                if (!is.na(as.numeric(input$rankmaxy))) 
                  maxy = as.numeric(input$rankmaxy)
                else maxy = sum(pctl)
                return(ruv_rankplot(fit, pctl, X.col = pvar) + 
                  coord_cartesian(xlim = c(0, maxx), ylim = c(0, 
                    maxy)))
            }))
            output$pplot = shiny::renderPlot(try({
                return(ggpplot() + recolor())
            }))
            output$varplot = shiny::renderPlot(try({
                if (is.null(input$bvar)) 
                  return(NULL)
                if (input$bvar == "") 
                  return(NULL)
                return(ruv_varianceplot(fit_summary(), X.col = input$bvar) + 
                  recolor())
            }))
            output$volcano = shiny::renderPlot(try({
                if (is.null(input$bvar)) 
                  return(NULL)
                if (input$bvar == "") 
                  return(NULL)
                return(ruv_volcano(fit_summary(), X.col = input$bvar) + 
                  recolor())
            }))
            output$ppbrush_info = shiny::renderTable({
                ppbrush = input$pplot_brush
                if (is.null(ppbrush)) 
                  return(NULL)
                x = ggpplot()$data$pplot.x
                y = ggpplot()$data$pplot.y
                keep = (x > ppbrush$xmin) & (x < ppbrush$xmax) & 
                  (y > ppbrush$ymin) & (y < ppbrush$ymax)
                if (sum(keep) == 0) 
                  return(NULL)
                selected = fit_summary()$C
                selected = cbind(data.frame(rownames = rownames(selected)), 
                  selected)
                selected = selected[keep, input$infocols]
                for (j in 1:ncol(selected)) if (!is.numeric(selected[, 
                  j]) & !is.logical(selected[, j])) 
                  selected[, j] = google_search(selected[, j])
                return(selected)
            }, sanitize.text.function = function(x) x)
            output$volbrush_info = shiny::renderTable({
                volbrush = input$volcano_brush
                if (is.null(volbrush)) 
                  return(NULL)
                b = fit_summary()$C[, paste0("b_", input$bvar)]
                p = fit_summary()$C[, paste0("p_", input$bvar)]
                volkeep = (b > volbrush$xmin) & (b < volbrush$xmax) & 
                  (p > 10^(-volbrush$ymax)) & (p < 10^(-volbrush$ymin))
                if (sum(volkeep) == 0) 
                  return(NULL)
                selected = fit_summary()$C
                selected = cbind(data.frame(rownames = rownames(selected)), 
                  selected)
                selected = selected[volkeep, c(paste0("b_", input$bvar), 
                  paste0("p_", input$bvar), input$infocols)]
                for (j in 1:ncol(selected)) if (!is.numeric(selected[, 
                  j]) & !is.logical(selected[, j])) 
                  selected[, j] = google_search(selected[, j])
                return(selected)
            }, sanitize.text.function = function(x) x)
            output$varbrush_info = shiny::renderTable({
                varbrush = input$varplot_brush
                if (is.null(varbrush)) 
                  return(NULL)
                b = fit_summary()$C[, paste0("b_", input$bvar)]
                s = fit_summary()$C[, "sigma2"]
                varkeep = (s > varbrush$xmin^4) & (s < varbrush$xmax^4) & 
                  (b^2 > varbrush$ymin^4) & (b^2 < varbrush$ymax^4)
                if (sum(varkeep) == 0) 
                  return(NULL)
                selected = fit_summary()$C
                selected = cbind(data.frame(rownames = rownames(selected)), 
                  selected)
                selected = selected[varkeep, c("sigma2", paste0("b_", 
                  input$bvar), input$infocols)]
                for (j in 1:ncol(selected)) if (!is.numeric(selected[, 
                  j]) & !is.logical(selected[, j])) 
                  selected[, j] = google_search(selected[, j])
                return(selected)
            }, sanitize.text.function = function(x) x)
            output$table = shiny::renderDataTable({
                if (is.null(fit_summary()) | is.null(input$infocols) | 
                  is.null(input$statcols)) 
                  return(NULL)
                toptable = fit_summary()$C
                toptable = cbind(data.frame(rownames = rownames(toptable)), 
                  toptable)
                if ("rownames" %in% input$infocols) 
                  keepcols = c("rownames", input$statcols, input$infocols[-1])
                else keepcols = c(input$statcols, input$infocols)
                toptable = toptable[, keepcols, drop = FALSE]
                for (j in 1:ncol(toptable)) if (!is.numeric(toptable[, 
                  j]) & !is.logical(toptable[, j])) 
                  toptable[, j] = google_search(toptable[, j])
                return(toptable)
            }, escape = FALSE)
            output$datamatrixtable = shiny::renderDataTable({
                Y()[, 1:min(ncol(Y()), 10)]
            })
            output$rowtable = shiny::renderDataTable({
                rowinfo()
            })
            output$coltable = shiny::renderDataTable({
                colinfo()
            })
            Y = shiny::reactive({
                if (is.null(input$Y.file)) 
                  return(Y.input)
                return(read.csv(input$Y.file$datapath, row.names = 1))
            })
            colinfo = shiny::reactive({
                if (is.null(input$colinfo.file)) 
                  return(colinfo.input)
                return(read.csv(input$colinfo.file$datapath, 
                  row.names = 1))
            })
            rowinfo = shiny::reactive({
                if (is.null(input$rowinfo.file)) 
                  return(rowinfo.input)
                return(read.csv(input$rowinfo.file$datapath, 
                  row.names = 1))
            })
            rowinfo3 = shiny::reactive({
                if (is.null(input$M) | !input$avg.III) 
                  return(rowinfo())
                if (length(input$burst) > 0) 
                  M = replicate.matrix(rowinfo()[, input$M], 
                    burst = input$burst)
                else M = replicate.matrix(rowinfo()[, input$M])
                return(collapse.replicates(rowinfo(), M))
            })
            Y1 = shiny::reactive({
                Y1 = Y()
                if (is.null(input$centercol) | is.null(input$centerrow)) 
                  return(NULL)
                if (input$centercol) 
                  Y1 = scale(Y1, scale = FALSE)
                if (is.null(input$ctl.I)) 
                  return(Y1)
                if (input$ctl.I == "All") 
                  ctl = rep(TRUE, ncol(Y1))
                else ctl = colinfo()[, which(names(colinfo()) == 
                  input$ctl.I)]
                if (input$centerrow) 
                  Y1 = RUV1(Y1, 1, ctl)
                return(Y1)
            })
            Y3 = shiny::reactive({
                Y3 = Y1()
                if (is.null(input$ctl.III)) 
                  return(NULL)
                if (input$ctl.III == "All") 
                  ctl = rep(TRUE, ncol(Y1()))
                else ctl = colinfo()[, which(names(colinfo()) == 
                  input$ctl.III)]
                if (is.null(input$M) | is.null(input$k.III.select) | 
                  is.null(input$avg.III)) 
                  return(NULL)
                if (length(input$M) == 0) 
                  return(NULL)
                if (length(input$burst) > 0) 
                  M = replicate.matrix(rowinfo()[, input$M], 
                    burst = input$burst)
                else M = replicate.matrix(rowinfo()[, input$M])
                if (input$k.III.select == "Maximum Possible") 
                  k = min(sum(ctl) - as.numeric(input$centerrow & 
                    input$ctl.I == input$ctl.III), nrow(Y3) - 
                    ncol(M))
                else k = input$k.III
                Y3 = RUVIII(Y3, M, ctl, k, average = input$avg.III)
                return(Y3)
            })
            Y1svd = shiny::reactive({
                if (is.null(Y1())) 
                  return(NULL)
                else return(svd(Y1()))
            })
            Y3svd = shiny::reactive({
                if (is.null(Y3())) 
                  return(NULL)
                else return(svd(Y3()))
            })
            output$burst = shiny::renderUI({
                if (length(input$M) == 0) 
                  return(NULL)
                shiny::checkboxGroupInput("burst", "Burst", choices = colnames(replicate.matrix(rowinfo()[, 
                  input$M])), selected = FALSE)
            })
            fit_summary = shiny::reactive({
                if (is.null(input$DE.Y.var)) 
                  return(NULL)
                if (input$DE.Y.var == "RUVIII") 
                  Y = Y3()
                else if (input$DE.Y.var == "Centered") 
                  Y = Y1()
                else Y = Y()
                rowinfo = rowinfo()
                if (!is.null(input$avg.III)) {
                  if (input$avg.III & input$DE.Y.var == "RUVIII") 
                    rowinfo = rowinfo3()
                }
                if (is.null(input$X) | is.null(input$ctl)) 
                  return(NULL)
                X = rowinfo[, which(names(rowinfo) == input$X)]
                if (input$ctl == "All") 
                  ctl = rep(TRUE, ncol(Y))
                else ctl = colinfo()[, which(names(colinfo()) == 
                  input$ctl)]
                if (!is.logical(ctl)) 
                  return(NULL)
                method = input$DE.method
                k = input$k
                if (is.null(k)) 
                  k = 0
                if (is.null(ctl) | is.null(X) | is.null(method)) 
                  return(NULL)
                if (method == "Unadjusted") 
                  fit = try(RUV2(Y, X, ctl, 0))
                if (method == "RUV2") 
                  fit = try(RUV2(Y, X, ctl, k))
                if (method == "RUV4") 
                  fit = try(RUV4(Y, X, ctl, k))
                if (method == "RUVinv") 
                  fit = try(RUVinv(Y, X, ctl))
                if (method == "RUVrinv") 
                  fit = try(RUVrinv(Y, X, ctl))
                ebayes = input$ebayes
                if (input$vtype == "Standard") 
                  var.type = "standard"
                if (input$vtype == "E. Bayes (Limma)") 
                  var.type = "ebayes"
                if (input$vtype == "Pooled") 
                  var.type = "pooled"
                if (input$ptype == "Standard") 
                  p.type = "standard"
                if (input$ptype == "Rescaled") 
                  p.type = "rsvar"
                if (input$ptype == "Empirical") 
                  p.type = "evar"
                if (class(fit) == "try-error") 
                  stop("Error encountered.  Too few negative controls and / or k too high?")
                fit = try(ruv_summary(Y, fit, rowinfo, colinfo(), 
                  var.type = var.type, p.type = p.type))
                if (class(fit) == "try-error") 
                  return(NULL)
                return(fit)
            })
            quick_fit = shiny::reactive({
                if (is.null(input$X)) 
                  return(NULL)
                X = rowinfo()[, which(names(rowinfo()) == input$X)]
                ctl = rep(TRUE, ncol(Y()))
                fit = RUV2(Y(), X, ctl, 0, inputcheck = FALSE)
                fit = ruv_summary(Y(), fit, rowinfo(), colinfo())
                return(fit)
            })
            recolor = shiny::reactive({
                if (is.null(fit_summary())) 
                  return(NULL)
                colorvar = input$color.by
                if (colorvar == "None") 
                  return(NULL)
                if (colorvar == "Negative Controls") 
                  colorvar = "fit.ctl"
                rval = list(aes_string(color = colorvar))
                if (!input$ccolor) 
                  return(rval)
                clevels = levels(as.factor(fit_summary()$C[, 
                  colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(rval)
                colorvals = rep(NA, nc)
                for (i in 1:nc) colorvals[i] = input[[paste0(make.names(clevels[i]), 
                  "_color")]]
                rval[[2]] = scale_color_manual(name = input$color.by, 
                  values = colorvals)
                return(rval)
            })
            recolor.row = shiny::reactive({
                if (is.null(input$color.by.row)) 
                  return(NULL)
                colorvar = input$color.by.row
                if (colorvar == "None") 
                  return(NULL)
                rval = list(aes_string(color = colorvar))
                if (!input$ccolor.row) 
                  return(rval)
                clevels = levels(as.factor(rowinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(rval)
                colorvals = rep(NA, nc)
                for (i in 1:nc) colorvals[i] = input[[paste0(make.names(clevels[i]), 
                  "_color.row")]]
                rval[[2]] = scale_color_manual(name = input$color.by.row, 
                  values = colorvals)
                return(rval)
            })
            recolor.rle = shiny::reactive({
                if (is.null(input$color.by.row)) 
                  return(NULL)
                colorvar = input$color.by.row
                if (colorvar == "None") 
                  return(aes_string(fill = NULL))
                rval = list(aes_string(fill = colorvar))
                if (!input$ccolor.row) 
                  return(rval)
                clevels = levels(as.factor(rowinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(rval)
                colorvals = rep(NA, nc)
                for (i in 1:nc) colorvals[i] = input[[paste0(make.names(clevels[i]), 
                  "_color.row")]]
                rval[[2]] = scale_fill_manual(name = input$color.by.row, 
                  values = colorvals)
                return(rval)
            })
            recolor.col = shiny::reactive({
                if (is.null(input$color.by.col)) 
                  return(NULL)
                colorvar = input$color.by.col
                if (colorvar == "None") 
                  return(NULL)
                rval = list(aes_string(color = colorvar))
                if (!input$ccolor.col) 
                  return(rval)
                clevels = levels(as.factor(colinfo()[, colorvar]))
                nc = length(clevels)
                if (nc > 12) 
                  return(rval)
                colorvals = rep(NA, nc)
                for (i in 1:nc) colorvals[i] = input[[paste0(make.names(clevels[i]), 
                  "_color.col")]]
                rval[[2]] = scale_color_manual(name = input$color.by.col, 
                  values = colorvals)
                return(rval)
            })
            svd.left = shiny::reactive({
                if (is.null(input$shape.row)) 
                  return(list(recolor.row(), aes_string(size = 4), 
                    scale_size_identity()))
                shapevar = input$shape.row
                if (shapevar == "None") 
                  return(list(recolor.row(), aes_string(size = 4), 
                    scale_size_identity()))
                return(list(recolor.row(), aes_string(shape = shapevar, 
                  size = 4), scale_size_identity()))
            })
            svd.left3 = shiny::reactive({
                if (!input$avg.III) 
                  return(svd.left())
                rval = list(aes_string(size = 4), scale_size_identity())
                if (input$color.by.row %in% colnames(rowinfo3())) 
                  rval = c(rval, list(recolor.row()))
                if (is.null(input$shape.row)) 
                  return(rval)
                shapevar = input$shape.row
                if (shapevar == "None") 
                  return(rval)
                if (input$shape.row %in% colnames(rowinfo3())) 
                  rval = c(rval, list(aes_string(shape = shapevar)))
                return(rval)
            })
            ggpplot = shiny::reactive({
                if (is.null(fit_summary()) | is.null(input$bvar)) 
                  return(NULL)
                if (input$bvar == "") 
                  return(NULL)
                return(ruv_projectionplot(fit_summary(), X.col = input$bvar, 
                  factor = input$ppfactor, adjusted = input$ppadjusted))
            })
            plot.k = shiny::reactive({
                return(as.numeric(input$plot.k))
            })
        }, options = options)
}
