prettify_percentage <- function(x, round_digits = 0L)
{
  round.digits <- as.integer(round_digits)
  return(sprintf("%+0.*f%%", round_digits, 100 * x))
}

prettify_single_number <- function(x, letter, round_digits)
{
  if (is.na(x) && !is.nan(x)) {
    return("NA")
  } else if (!is.finite(x)) {
    return(as.character(x))
  } else if ((letter == "" && abs(x) >= 1e9) || letter == "B") {
    return(sprintf("%0.*fB", round_digits, x / 1e9))
  } else if ((letter == "" && abs(x) >= 1e6) || letter == "M") {
    return(sprintf("%0.*fM", round_digits, x / 1e6))
  } else if ((letter == "" && abs(x) >= 1e3) || letter == "K") {
    return(sprintf("%0.*fK", round_digits, x / 1e3))
  } else if (abs(x) >= 1 || x == 0) {
    return(sprintf("%0.*f", round_digits, x))
  } else {
    # Calculate position of first non-zero digit after the decimal point
    first.nonzero <- -floor(log10(abs(x)))
    return(sprintf("%0.*f", round_digits + first.nonzero - 1, x))
  }
}

prettify_number <- function(x, letter = "", round_digits = 1L)
{
  round_digits <- as.integer(round_digits[1])
  letter <- rep(letter, length.out = length(x))
  output <- sapply(seq_along(x), function(index) {
    prettify_single_number(x[index], letter[index], round_digits)
    })
  return(output)
}

string_trim <- function(x)
{ 
  gsub("^\\s+|\\s+$", "", x)
}

format_number <- function(x, digits)
{
  string_trim(format(x, digits = digits))
}

format_percent <- function(x, digits) 
{
  string_trim(paste0(format(x * 100, digits = digits), "%"))
}

format_confidence_intervals <- function(a, b, digits)
{
  paste0("[", string_trim(format(a, digits = min(digits, 2))),
         ", ", string_trim(format(b, digits = min(digits, 2))),
         "]")
}

format_percent_confidence_intervals <- function(a, b, digits)
{
  paste0("[", string_trim(format(a * 100, digits = min(digits, 2))),
         "%, ", string_trim(format(b * 100, digits = min(digits, 2))),
         "%]")
}

identify_number_abbreviation <- function(abbreviated_number)
{
  letter <- substr(abbreviated_number, nchar(abbreviated_number),
                   nchar(abbreviated_number))
  letter[!is.element(letter, c("B", "M", "K"))] <- "none"
  return(letter)
}
