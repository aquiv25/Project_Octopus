library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

grDevices::pdf.options(family = "Helvetica")

# ── Load data ──────────────────────────────────────────────────────────────────
df_raw <- read.delim(
  "/Users/dassagaripova/Downloads/enrichment.all.tsv",
  comment.char = "#", header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, quote = ""
)
colnames(df_raw) <- c("category","term_id","term","obs","bg","strength","signal","fdr",
                       "protein_ids","protein_labels")

df_raw$fdr        <- as.numeric(df_raw$fdr)
df_raw$obs        <- as.numeric(df_raw$obs)
df_raw$bg         <- as.numeric(df_raw$bg)
df_raw$log10fdr   <- -log10(df_raw$fdr)
df_raw$gene_ratio <- df_raw$obs / df_raw$bg * 100

# ── Top 7 per category ────────────────────────────────────────────────────────
TOP_N     <- 7
cat_order_full <- c("GO Process","GO Component","GO Function",
                    "KEGG","Reactome","InterPro",
                    "UniProt Keywords","STRING clusters","SMART","COMPARTMENTS")
available   <- unique(df_raw$category)
cat_order   <- cat_order_full[cat_order_full %in% available]
cat(sprintf("Categories: %s\n", paste(cat_order, collapse = ", ")))

df <- df_raw %>%
  filter(category %in% cat_order) %>%
  group_by(category) %>%
  slice_min(order_by = fdr, n = TOP_N) %>%
  ungroup()

df$category <- factor(df$category, levels = rev(cat_order))

# ── Wrap long labels ──────────────────────────────────────────────────────────
wrap_label <- function(x, width = 45) {
  sapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"))
}
df$term_wrapped <- wrap_label(df$term, 45)

df <- df %>% arrange(category, log10fdr)
df$term_wrapped <- factor(df$term_wrapped, levels = unique(df$term_wrapped))

# ── Color palette — Z-DNA / red/crimson ───────────────────────────────────────
cat_colors_all <- c(
  "GO Process"         = "#B71C1C",
  "GO Component"       = "#C62828",
  "GO Function"        = "#D32F2F",
  "KEGG"               = "#880E0E",
  "Reactome"           = "#AD1457",
  "InterPro"           = "#6A0000",
  "UniProt Keywords"   = "#E53935",
  "STRING clusters"    = "#4A0000",
  "SMART"              = "#FF5252",
  "COMPARTMENTS"       = "#FF8A80"
)
cat_colors <- cat_colors_all[cat_order]

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = log10fdr, y = term_wrapped,
                    size = gene_ratio, fill = category)) +

  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = category),
            alpha = 0.06, inherit.aes = FALSE, show.legend = FALSE) +

  geom_vline(xintercept = -log10(0.05), color = "#888888",
             linetype = "dashed", linewidth = 0.4) +

  geom_point(shape = 21, color = "white", stroke = 0.25, alpha = 0.90) +

  geom_text(data = df %>% group_by(category) %>% slice_min(fdr, n = 1),
            aes(label = formatC(log10fdr, digits = 1, format = "f")),
            size = 2.2, color = "white", fontface = "bold") +

  scale_fill_manual(values = cat_colors, name = "Category",
                    guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_size_continuous(
    name   = "Gene ratio (%)",
    range  = c(1.5, 12),
    breaks = c(1, 5, 10, 20),
    labels = c("1%", "5%", "10%", "20%")
  ) +
  scale_x_continuous(
    name   = expression(bold(-log[10](FDR))),
    expand = expansion(mult = c(0.01, 0.07)),
    labels = label_number(accuracy = 1)
  ) +

  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +

  theme_minimal(base_size = 10) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "#FFF8F8", color = NA),
    panel.grid.major.x = element_line(color = "#FFCDD2", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.4, "lines"),
    panel.border       = element_rect(color = "#FFCDD2", fill = NA, linewidth = 0.4),

    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                      size = 9, face = "bold",
                                      margin = margin(r = 8)),
    strip.background   = element_blank(),

    axis.title.x  = element_text(size = 10, face = "bold",
                                  color = "#7F0000", margin = margin(t = 8)),
    axis.title.y  = element_blank(),
    axis.text.x   = element_text(size = 9, color = "#333333"),
    axis.text.y   = element_text(size = 8.5, color = "#111111",
                                  hjust = 1, lineheight = 0.9),
    axis.ticks.x  = element_line(color = "#AAAAAA", linewidth = 0.4),
    axis.ticks.y  = element_blank(),

    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.title      = element_text(size = 9, face = "bold", color = "#7F0000"),
    legend.text       = element_text(size = 8.5),
    legend.key.size   = unit(1.1, "lines"),
    legend.margin     = margin(t = 8),
    legend.spacing.x  = unit(0.4, "cm"),

    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5,
                                  color = "#7F0000", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "#555555",
                                  margin = margin(b = 10)),
    plot.caption  = element_text(size = 7.5, color = "#777777",
                                  hjust = 1, margin = margin(t = 8)),
    plot.margin   = margin(12, 20, 10, 8)
  ) +

  labs(
    title    = "Functional Enrichment — Z-DNA (Top 7 terms per category)",
    subtitle = "Bubble size = gene ratio (obs/bg)  ·  Dashed line: FDR = 0.05",
    caption  = "Source: STRING enrichment · Octopus bimaculoides · Z-DNABERT & Z-Hunter predictions"
  )

# ── Save ──────────────────────────────────────────────────────────────────────
n_terms <- nrow(df)
n_lines <- sum(sapply(as.character(df$term_wrapped),
                       function(x) length(strsplit(x,"\n")[[1]])))
plot_h  <- max(8, n_lines * 0.22 + n_terms * 0.1 + 4)
plot_w  <- 14

out_png <- "/Users/dassagaripova/Downloads/project/enrichment_bubble_zdna.png"
out_svg <- "/Users/dassagaripova/Downloads/project/enrichment_bubble_zdna.svg"

ggsave(out_png, p, width = plot_w, height = plot_h, dpi = 300)
tryCatch(ggsave(out_svg, p, width = plot_w, height = plot_h, device = "svg"),
         error = function(e) message("SVG skipped"))

cat(sprintf("✅ Z-DNA PNG: %s  (%.1f × %.1f in, %d terms)\n",
            out_png, plot_w, plot_h, n_terms))
