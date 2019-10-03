library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
library(broom)
library(gridExtra)

source("c:\\r\\intensity\\INT\\Intensity_functions.R") # load functions stored in the external file

path <- "c:\\r\\intensity\\data4\\" # directory containing the data Excel files

strainCodes <- c("A", "B", "C", "E", "F", "J", "M", "P", "R", "43") # list of strain codenames included in the experiment

strainNames <- c(  # list of actual gene names
  "RSP4/6",
  "RSP2",
  "RSP3",
  "IAD",
  "OAD2",
  "FAP20",
  "DRC4",
  "PFR2",
  "Rib72",
  "CFAP43"
)


fileNames <- paste0(path, strainCodes, ".xlsx") # complete names and paths of the files

hours <- c(0, 2, 4, 16, 24, 48) # list of timepoints of the experiment

pixelSize <- 0.065

fileNames %>% # take the sequence of the filenames

  set_names(strainCodes) %>% # assign them the strain codenames

  map_df(~ readFile(.), .id = "Line") %>% # run  two nested functions (readFile, cleaner_sheet) to get a table of raw data

  mutate(Hour = as.numeric(substr(Hour, 1, nchar(Hour) - 1))) %>% # change the hour format to numeric

  group_by(Line, Hour, Cell) %>% # group by individual cells``

  filter(n() > 101) %>% # filter away very short ones

  {
    . ->> rawData
  } # %>% # keep the intermediate - raw data

## ---------------- proces raw data - find decaypoint, fit models etc. ------------

rawData %>%
  mutate(
    Point = findPoints(Red, Green), # [1]~1, [2]~beginning of flagellum, [3]~ decaypoint,[4]~devoidpoint, [5]end of flagellum, [6]end, rest are 0's
    LineName = strainNames[which(strainCodes == first(Line))], # assign the actual gene names
    RedIntNorm = Red / max(Red), # add normalized red intensity (0-1)
    GreenIntNorm = Green / max(Green), # add normalized green intensity (0-1)
    TotalLength = Position[Point[5] ] - Position [Point[2] ],
    DecayLength = Position[Point[5]] - Position [Point[3]],
    DevoidLength = Position[Point[5]] - Position [Point[4]],
    Area = sum(GreenIntNorm[Point[3]:Point[5]] -
      RedIntNorm[Point[3]:Point[5]]),
    FitAlpha = as.numeric(tryTidyFit(RedIntNorm, Point[3])), # tryTidyFit function returns vector, where 1st elements is the alpha, rest is the prediction
    Alpha = FitAlpha[1], # separate alpha
    Fit = c(NA, FitAlpha[-1]), # and the predictions
    MeanIntensity = mean(Red[round((max(1, Point[3] - (5 / pixelSize))):max(3 / pixelSize, Point[3] - (2 / pixelSize)))])
  ) %>%
  select(-FitAlpha) %>% # remove the temporary column

  unite("Cell_ID", c(Line, Hour, Cell), remove = F) %>% # create unique identifier for each cell

  filter(!Cell_ID %in% c("43_48_cell19", "A_48_cell24", "P_16_cell07")) %>% # fishy cells

  ungroup() %>%
  {
    . ->> fullData
  } %>% # keep the intermediate - full data

  nest(c(Red, Green, Position, Point, Fit, RedIntNorm, GreenIntNorm)) %>% # pack the longer data, so that each cell contained on one row

  {
    . ->> fullDataNested
  }

## ------------- Plot all the individual curves, this takes more than 7 minutes ---------------

fullDataNested %>%
  unnest() %>%
  group_by(Cell_ID) %>%
  group_map(~ {
    pl <- ggplot(.x, aes(Position, RedIntNorm)) + # Prepares a plot for each group, returns all of them as a list
      geom_line(colour = "red", lwd = 1.2) +
      geom_line(
        data = .x[round((max(1, .x$Point[3] - (5 / pixelSize)))):
        round(max(3 / pixelSize, .x$Point[3] - (2 / pixelSize))), ],
        aes(
          x = Position,
          y = RedIntNorm
        ),
        lwd = 1,
        colour = "blue"
      ) +
      geom_line(aes(y = Fit), lwd = 1.2, na.rm = T) +
      geom_vline(xintercept = .x$Point[1:6] * pixelSize) +
      geom_vline(xintercept = .x$Point[4] * pixelSize, colour = "green") +
      geom_vline(xintercept = .x$Point[3] * pixelSize, colour = "blue") +
      geom_line(aes(Position, GreenIntNorm), colour = "darkgreen", alpha = 0.6) +
      labs(title = paste(
        "Line", .x$Line, "Hour", .x$Hour, .x$Cell,
        "Alpha", round(.x$Alpha, 2),
        "Mean Intensity", round(.x$MeanIntensity, 2),
        "Decay length", round(.x$DecayLength, 2)
      ))
  }) -> listOfPlots

ggsave("c:\\r\\intensity\\data4\\allPlots.pdf",
  marrangeGrob(grobs = listOfPlots, nrow = 2, ncol = 1),
  width = 20, height = 20, units = "cm"
)
