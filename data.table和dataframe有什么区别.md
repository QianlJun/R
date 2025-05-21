#  data.table和dataframe有什么区别

data.table按照列名称选取数据时，要在变量前加".."

举例：

```
 # 选取需要的列
 cols_to_select <- c("variant", case_ac_col, case_an_col, ctr_ac_col, ctr_an_col)
 variant_data <- variant_data[, ..cols_to_select]
```

