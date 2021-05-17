# from ?lm (reformatted a bit)
d_d9 <-
  data.frame(
    weight=c(c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14),
             c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)),
    group=gl(2, 10, 20, labels = c("Ctl","Trt"))
  )
test_lm_model <- lm(weight ~ group, data=d_d9)

test_that("add_predictions2", {
  expect_equal(
    add_predictions2(data=d_d9, model=test_lm_model)$pred,
    unname(predict(test_lm_model))
  )
  expect_equal(
    add_predictions2(data=d_d9, model=test_lm_model, var="newpred")$newpred,
    unname(predict(test_lm_model)),
    info="var is respected"
  )
  expect_s3_class(
    add_predictions2(data=d_d9, model=test_lm_model, var="newpred"),
    "data.frame"
  )
})
