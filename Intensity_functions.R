#------- find all local maxima ----------
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#------- find all local minima ----------
localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#------- find decay point using the red and green difference method ----------
bpDiff <- function(sampleGreenNormMax, sampleRedNormMax) {
  maxDiff <- sampleGreenNormMax - sampleRedNormMax
  maxDiff.plot <- maxDiff
  maxDiff[maxDiff >= 0] <- 0
  maxDiff[maxDiff < 0] <- 1
  maxDiffrev <- rev(abs(maxDiff)) # reverse the order of the diffrence series
  chpt <- max(which((cumsum(maxDiffrev) < 5))) # changepoint of the difference series is a point where the series is at least 5 times below0
  breakpoint <- length(maxDiff) - chpt + 1

  if (maxDiff.plot[breakpoint] < 0) {
    if (any(maxDiff.plot[breakpoint:length(maxDiff.plot)] >= 0)) {
      breakpoint <- min(which(maxDiff.plot[breakpoint:length(maxDiff.plot)] >= 0) + breakpoint - 1)
    }

    else {
      breakpoint <- length(maxDiff.plot)
    }
  }

  return(breakpoint)
}

#------- find decay point using the cumulative sum method  ----------
bpCusum <- function(sampleGreenNormSum, sampleRedNormSum) {
  cSRed <- cumsum(sampleRedNormSum)
  cSGreen <- cumsum(sampleGreenNormSum)
  # red signal + local differences between red and green cumulative cum curve
  rozdil <- (cSRed + (cSRed - cSGreen))
  # flattened rozdil curve
  corrected.rozdil <- rozdil - (seq_along(rozdil) / max(seq_along(rozdil)))

  if (max(corrected.rozdil) >= 0.3) breakpoint <- which(corrected.rozdil == max(corrected.rozdil))

  if (max(corrected.rozdil) < 0.3) breakpoint <- max(localMaxima(corrected.rozdil))

  return(breakpoint)
}



#------ try to fit assymptote to the end of the red signal ----------
tryFit <- function() {
  fit <- findFit(breakpoint_find, EOF, sampleRedNormMax)
  redTailX <- breakpoint_find:(EOF + 50)
  am <- tidy(fit)
  alpha_exp <- am[3, 2]
  if (is.na(alpha_exp)) alpha_exp <- -10 # if all fails, hard set the alpha to -10, assuming linear gradient
}



# -----non-linear regression to quantify the red signal decay-----
findFit <- function(breakpoint, EOF, sampleRedNormMax) {
  redTailX <- breakpoint:(EOF + 50)

  # print(length(redTailX))

  redTailY <- c(sampleRedNormMax[breakpoint:EOF], rep(0, 50))

  df <- data.frame(x = redTailX, y = redTailY)

  fit <- tryCatch({
    fit <- nls(y ~ SSasymp(x, yf, y0, log_alpha), data = df)
  },
  error = function(e) {
    fit <- lm(y ~ x, data = df)
    return(fit)
  },

  warning = function(w) {
    fit <- lm(y ~ x, data = df)
    return(fit)
  }
  )

  return(fit)
}


## function for group sampling
## downloaded from https://cmdlinetips.com/2019/07/how-to-randomly-select-groups-in-r-with-dplyr/

sample_n_groups <- function(grouped_df, size, replace = FALSE, weight = NULL) {
  grp_var <- grouped_df %>%
    groups() %>%
    unlist() %>%
    as.character()
  random_grp <- grouped_df %>%
    summarise() %>%
    sample_n(size, replace, weight) %>%
    mutate(unique_id = 1:NROW(.))
  grouped_df %>%
    right_join(random_grp, by = grp_var) %>%
    group_by_(grp_var)
}

write_sheets_as_csv <- function(sheet, path) {
  pathbase <- path %>%
    tools::file_path_sans_ext()
  path %>%
    read_excel(sheet = sheet) %>%
    write_csv(paste0(pathbase, "-", sheet, ".csv"))
}




#---- select just the needed columns out of an Excel sheet-------
cleaner_sheet <- function(sheet) {
  listOfcellsIndex <- which(grepl("cell", sheet[2, ])) # finds columns containing cell identifier, which is expected to be in the second row

  # produces long table of fluorescence intensity observations
  # needed columns are localized by their relative position to the "cellxx" identifier
  rst <- tibble(
    Cell = unlist(rep(sheet[2, listOfcellsIndex], each = (nrow(sheet) - 3))),
    Position = as.numeric(unlist(sheet[4:nrow(sheet), listOfcellsIndex + 4])),
    Red = as.numeric(unlist(sheet[4:nrow(sheet), listOfcellsIndex + 2])),
    Green = as.numeric(unlist(sheet[4:nrow(sheet), listOfcellsIndex + 6]))
  ) %>%
    drop_na() # shorter columns are filled with NAs to match the longest column, this removes those lines
}

#---- read sheets in an  Excel file-------
readFile <- function(fileName) {
  fileName %>%

    excel_sheets() %>% # get sheet names from the files

    keep(~ grepl(paste(hours, collapse = "|"), .x)) %>% # keep only sheets containing one of the hour numbers in their title

    set_names(str_extract_all(., "[0-9]+h")) %>% #  "[0-9]+" extract the hours from the sheet titles and assign as names

    map(~ read_excel(path = fileName, sheet = .x)) %>% # return list of excel sheets as tables

    map_df(~ cleaner_sheet(.), .id = "Hour") # calls the sheet cleaning function and binds the clean sheets to one long dataframe
}


#------------find points in the intensity profile (beginning, decay, end etc...)----------------
findPoints <- function(rr, gg) {
  sampleGreen <- gg
  sampleRed <- rr


  # normalize to max intensity = 1

  sampleGreenNormMax <- sampleGreen / max(sampleGreen)
  sampleRedNormMax <- sampleRed / max(sampleRed)


  # use drop of green signal to find the end of the flagellum (EOF)

  lowRatio <- (which((sampleGreen / mean(sampleGreen)) < 0.1)) # look for points where the green signal intensity drops below 0.1 of the mean

  if (length(lowRatio) > 0) {
    lowRatio <- lowRatio[lowRatio > 50] # we are not interested in values too close to the flagellum base

    EOF <- min(lowRatio, length(sampleGreen))
  } # out of the remaining, take the point closest to the base

  else {
    EOF <- length(sampleGreen) - 1
  } # if no alternative found, end of flagellum remains equal to the end of the data


  # find beginning of flagellum (BOF) as a value closest to the 25th percentile the first 100 points
  BOF <- min(which(
    sampleGreenNormMax [1:100] == quantile(sampleGreenNormMax[1:100],
      probs = 0.11,
      names = FALSE,
      type = 3
    ) # type 3 returns one of the actual values, not an interpolation
  ))

  if (is.infinite(BOF) | BOF == 1) BOF <- 2

  if (EOF == length(sampleGreen)) EOF <- length(sampleGreen) - 1

  #--------find breakpoints-----------

  breakpoint <- EOF
  breakpoint <- bpDiff(
    sampleGreenNormMax[BOF:EOF],
    sampleRedNormMax[BOF:EOF]
  ) # find breakpoint using the red/green difference method


  if ((EOF - breakpoint) < 90) { # choose which breakpoint method to use, bpC works better for longer signal decays

    # normalize to the total sum of intensity = 1

    sampleGreenNormSum <- sampleGreen / sum(sampleGreen)
    sampleRedNormSum <- sampleRed / sum(sampleRed)

    breakpoint <- bpCusum(
      sampleGreenNormSum [BOF:EOF],
      sampleRedNormSum[BOF:EOF]
    )
  }

  breakpoint <- breakpoint + BOF - 1

  redTailDensity <- cumsum(rev(sampleRed[breakpoint:EOF])) / seq_along(sampleRed[breakpoint:EOF])
  devoidPoint <- max(which(redTailDensity < 100), 1)
  devoidPoint <- EOF - devoidPoint + 1

  ptsVec <- c(1, BOF, breakpoint, devoidPoint, EOF, length(sampleGreen), rep(0, length(sampleGreen) - 6))

  return(ptsVec)
}


#---------- non-linear regression to quantify the red signal decay-------------

findTidyFit <- function(redTail) {
  redTailY <- c(redTail, rep(0, 50))

  df <- data.frame(x = seq_along(redTailY), y = redTailY)

  fit <- tryCatch({
    fit <- nls(y ~ SSasymp(x, yf, y0, log_alpha), data = df)
  },
  error = function(e) {
    fit <- lm(y ~ x, data = df)
    return(fit)
  },

  warning = function(w) {
    fit <- lm(y ~ x, data = df)
    return(fit)
  }
  )

  return(fit)
}

tryTidyFit <- function(Red, bP) {
  #------ try to fit assymptote to the end of the red signal ----------


  fit <- findTidyFit(Red[bP:length(Red)])

  am <- tidy(fit)
  alpha_exp <- am[3, 2]
  if (is.na(alpha_exp)) { # the first try was nor succesful, try again with longer Red sample
    locMax <- max(localMaxima(Red[1:bP])) # find local maximum of red closest to the decaypoint. . .
    fit <- findTidyFit(Red [ locMax:length(Red) ]) # . . . and fit from there
    am <- tidy(fit)
    alpha_exp <- am[3, 2]

    if (is.na(alpha_exp)) {
      return(c(-10, rep(NA, length(Red) - 1)))
    } # second attempt failed,  set the alpha to -10, assuming linear gradient
    else {
      pFit <- predict(fit) # second attempt succeful

      return(c(alpha_exp, rep(NA, length(Red) - (length(pFit[1:(length(pFit) - 50) ]) + 1)), pFit[1:(length(pFit) - 50) ]))
    }
  }
  else {
    pFit <- predict(fit) # first try was succesful
    return(c(alpha_exp, rep(NA, length(Red) - (length(pFit[1:(length(pFit) - 50) ]) + 1)), pFit[1:(length(pFit) - 50) ]))
  }
}
