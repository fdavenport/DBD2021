library(ggplot2);

## figures can be up to 7" or 7.2" wide for a double column
pub_font_size <- 7
extdata_title_size <- 9


## version of theme_classic where all dark gray colors have been replaced with true black
theme_classic_black <- function(base_size = 11, base_family = "",
                                base_line_size = base_size / 22,
                                base_rect_size = base_size / 22) {
    ## default line thickness used for all themes
    half_line <- base_size / 2

    theme_classic(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
    ) %+replace%
        theme(axis.text = element_text(colour = "black", size = rel(0.8)),
              axis.ticks = element_line(colour = "black", size = rel(0.5)),
              strip.text = element_text(colour = "black",
                                        size = rel(0.8),
                                        margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)),

              strip.background = element_blank(),

              complete = TRUE
    )
}

## map theme for publication-size figures
map_theme_pub <- theme_void() +
    theme(text = element_text(size = pub_font_size),
          plot.title = element_text(size = pub_font_size, hjust = 0.5),
          legend.text = element_text(size = pub_font_size),
          legend.title = element_text(size = pub_font_size),
          strip.background = element_blank(),
          strip.text = element_text(size = pub_font_size, hjust = 0.5),
          plot.margin = margin(t = 0, r = -2, b = -6, l = -2, unit = "pt"),
          legend.box.margin = margin(t = -6, r = 2, b = 5, l = -1, unit = "pt"))


## theme for publication-size figures
theme_pub <- theme_classic_black() +
    theme(text = element_text(size = pub_font_size),
          legend.text = element_text(size = pub_font_size),
          legend.title = element_text(size = pub_font_size),
          axis.text = element_text(size = pub_font_size),
          plot.title = element_text(size = pub_font_size, hjust = 0.5),
          strip.text = element_text(size = pub_font_size, hjust = 0.5))

## theme for extended data figures
theme_SI <- theme_classic_black() +
    theme(text = element_text(size = pub_font_size),
          legend.text = element_text(size = pub_font_size),
          legend.title = element_text(size = pub_font_size),
          axis.text = element_text(size = pub_font_size),
          plot.title = element_text(size = extdata_title_size, hjust = 0),
          strip.text = element_text(size = 8, hjust = 0.5),
          legend.key.size = unit(10, "pt"))

## ADD MAP FORMATTING
add_map_formatting <- function(plot_object, scale_name = waiver(), lab = waiver(),
                           val = brewer.pal(10, "BrBG"), drp = TRUE,
                           showlim = FALSE){
    plot_object +
        scale_fill_manual(name = scale_name,
                          values = val,
                          labels = lab,
                          drop = drp,
                          guide = guide_coloursteps(barwidth = 7,
                                                barheight = 0.3,
                                                ticks = TRUE,
                                                ticks.colour = "gray36",
                                                title.position = "top",
                                                title.hjust = 0.5,
                                                frame.colour = "gray36",
                                                frame.linewidth = .1,
                                                show.limits=showlim)) +
        coord_sf(crs = st_crs(2163), xlim = c(-2300000, 2500000),
                 ylim = c(-2100000,730000)) +
        map_theme_pub +
        theme(legend.position = "bottom",
              legend.title = element_text(size = pub_font_size),
              legend.text = element_text(size = pub_font_size),
              plot.title = element_text(size = pub_font_size),
              panel.grid.major=element_line(colour="transparent"),
              panel.spacing.y = unit(-3, "pt"),
              strip.text = element_text(angle = 270, size = pub_font_size,
                                        colour = "transparent"))
}

## -----------------------------------------------------------------------------
## OTHER PLOTTING VARIABLES
## -----------------------------------------------------------------------------

q <- c(0.5, 0.75, 0.95, 0.99)
quant_names <- paste0("q", q)
quant_labels <- paste0(str_pad(gsub("0.", "", q), width = 2,
                               pad = "0", side = "right"), "th percentile")
names(quant_labels) <- quant_names
quant_labeller <- as_labeller(quant_labels)

## CMIP forcing colors
forcing_colors <- c("nat" = "dodgerblue3",
                    "hist+\nRCP8.5" = "black",
                    "RCP2.6     " = rgb(196, 121, 0, maxColorValue = 255),
                    "  RCP8.5\n" = rgb(153, 0, 2, maxColorValue=255))

## REGION labels and colors
region_labels <- c("NW", "W", "SW", "NR", "S", "UM", "C", "SE", "NE")
region_label_long<- c("North\nwest", "West", "South\nwest", "Northern\nRockies\nPlains",
                      "South", "Upper\nMidwest", "Central", "Southeast", "Northeast")
region_spacing <- c(2, 4.5, 7.5, 12, 17.5, 22.5, 28, 34.5, 43)
region_colors <- c(brewer.pal(8, "Dark2"), "dodgerblue4")
names(region_colors) <- region_order[1:9]

season_col <- c("#33A02C", "#4292C6", "#FF7F00", "gray38")
