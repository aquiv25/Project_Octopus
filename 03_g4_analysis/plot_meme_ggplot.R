library(ggplot2)
library(ggseqlogo)

# Parse motif 1 matrix from meme.txt
txt <- readLines("meme_output/meme.txt")

# Find letter-probability matrix for motif 1
start <- grep("letter-probability matrix", txt)[1] + 1
rows  <- c()
i <- start
while (i <= length(txt)) {
  line <- trimws(txt[i])
  vals <- suppressWarnings(as.numeric(strsplit(line, "\\s+")[[1]]))
  if (length(vals) == 4 && !any(is.na(vals))) {
    rows <- rbind(rows, vals)
    i <- i + 1
  } else break
}

# Matrix: rows = positions, cols = A C G T
# ggseqlogo needs matrix: rows = letters (A C G T), cols = positions
pwm <- t(rows)
rownames(pwm) <- c("A", "C", "G", "T")
colnames(pwm) <- seq_len(ncol(pwm))

# Custom DNA colors matching our palette
my_colors <- c(A = "#2CA02C", C = "#1565C0", G = "#FFB300", T = "#D62728")

p <- ggseqlogo(pwm,
               method   = "bits",
               col_scheme = make_col_scheme(
                 chars  = c("A","C","G","T"),
                 cols   = unname(my_colors)
               )) +
  labs(
    title   = NULL,
    subtitle= NULL,
    x       = "Position",
    y       = "Information content (bits)",
    caption = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5,  hjust = 0.5, color = "#444"),
    plot.caption  = element_text(size = 9, color = "#777", hjust = 0.5),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    panel.grid.major.y = element_line(color = "#eeeeee"),
    plot.margin   = margin(12, 16, 12, 16)
  )

ggsave("meme_output/motif1_ggplot.png", p,
       width = 9, height = 4.5, dpi = 180, bg = "white")
cat("✅ meme_output/motif1_ggplot.png\n")
