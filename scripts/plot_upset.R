
        if (!require("UpSetR", quietly = TRUE)) install.packages("UpSetR", repos="http://cran.us.r-project.org")
        library(UpSetR)

        df <- read.csv("presence_absence_matrix.csv", check.names=FALSE)
        plot_sets <- c("LSDV", "CAPV", "GTPV", "ancient", "pre_modern")

        # Publication-ready Okabe-Ito inspired palette
        base_colors <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#999999")
        set_colors <- base_colors[1:length(plot_sets)]

        # Bumped dimensions (12x8) for large fonts
        png("upset_plot.png", width=12, height=8, units="in", res=300)
        p <- upset(df, sets=plot_sets, keep.order=TRUE, order.by="freq",
                   mainbar.y.label = "Number of Shared CDS",
                   sets.x.label = "Total CDS per Genome",
                   sets.bar.color = set_colors,
                   point.size = 4, line.size = 1.5,
                   text.scale = c(2.0, 1.8, 2.0, 1.5, 2.0, 1.6))
        print(p)
        dev.off()
        