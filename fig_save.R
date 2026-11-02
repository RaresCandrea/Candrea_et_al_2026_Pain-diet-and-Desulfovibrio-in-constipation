# Common figure export: TIFF 600 dpi (LZW) for submission, PNG 300 dpi for review (ragg), vector PDF (cairo).
# Note: italic taxon names in plotmath must be quoted strings, italic("Name"), not bare symbols.
save_fig <- function(fig, stem, W_mm = 174, H_mm) {
  W <- W_mm / 25.4; H <- H_mm / 25.4
  ggplot2::ggsave(paste0(stem, ".tiff"), fig, width = W, height = H, dpi = 600, device = ragg::agg_tiff,
                  compression = "lzw", bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), fig, width = W, height = H, dpi = 300, device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), fig, width = W, height = H, device = grDevices::cairo_pdf)
  # editable PowerPoint (one slide, figure at its final size; text, axes and marks are editable shapes)
  if (isTRUE(getOption("fig.pptx", FALSE))) {
    doc <- officer::read_pptx()
    doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
    doc <- officer::ph_with(doc, rvg::dml(ggobj = fig, fonts = list(sans = "Arial")),
                            location = officer::ph_location(left = 0.2, top = 0.2, width = W, height = H))
    print(doc, target = paste0(stem, ".pptx"))
    # slide size = figure size + 0.2 inch margin on each side (officer version without set_slide_size)
    system2("python3", c("-c", shQuote(sprintf(paste0(
      "from pptx import Presentation; from pptx.util import Inches; p=Presentation('%s.pptx'); ",
      "p.slide_width=Inches(%f); p.slide_height=Inches(%f); p.save('%s.pptx')"), stem, W + 0.4, H + 0.4, stem))))
  }
  invisible(stem)
}
