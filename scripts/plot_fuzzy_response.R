library(ggplot2)
library(plotly)

# ---- 1. 发射概率 ~ 杂合率 ----
H <- seq(0, 0.5, by = 0.01)
predicted_homogeneity <- 0.85 - 0.2 * (H / 0.5)
predicted_homogeneity <- pmax(0.7, pmin(0.95, predicted_homogeneity))

df1 <- data.frame(H, predicted_homogeneity)

p1 <- ggplot(df1, aes(x = H, y = predicted_homogeneity)) +
  geom_line(color = "#0073C2FF", size = 1.2) +
  geom_point(color = "#0073C2FF") +
  labs(
    title = "模糊控制下发射概率随杂合率变化",
    x = "杂合率 H",
    y = "predicted_homogeneity"
  ) +
  theme_minimal(base_size = 14)

# ---- 2. 转移概率基准 ~ SNP密度与样本数 ----
# ---- 2. 转移概率基准 ~ SNP密度与样本数 ----
D <- seq(0, 1000, length.out = 40)  # SNP/Mb
N <- seq(50, 1000, length.out = 40) # 样本数
grid <- expand.grid(D = D, N = N)

# 模糊规则：高D、高N => cross_count高
base_cross <- 300 +
  300 * pmin(pmax((grid$D - 300)/400, 0), 1) +  # D效应
  200 * pmin(pmax((grid$N - 400)/600, 0), 1)    # N效应
grid$predicted_cross_count <- pmax(100, pmin(1200, base_cross))

# 将长表转为矩阵格式（行 = N, 列 = D）
z_matrix <- matrix(grid$predicted_cross_count,
                   nrow = length(N), ncol = length(D), byrow = TRUE)

# 绘制 3D 曲面
p2 <- plot_ly(
  x = D, y = N, z = z_matrix,
  type = "surface",
  colorscale = list(c(0, 1), c("skyblue", "orange")),
  contours = list(z = list(show = TRUE))
) %>%
  layout(
    title = "模糊控制下 HMM 断点数响应曲面",
    scene = list(
      xaxis = list(title = "SNP 密度 D (个/Mb)"),
      yaxis = list(title = "样本数 N"),
      zaxis = list(title = "predicted_cross_count")
    )
  )

# ---- 3. 显示图 ----
print(p1)
p2
