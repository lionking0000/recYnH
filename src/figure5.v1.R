# Assessment of the floral origin of honey via proteomic tools


layer1=c(0,1.34,0,0.56,0.23,1,2.7,0.6,3.57,0.34,0.71,0.28,0.51,1.71,0.79,0.07,1.47,0.03,2.77,0.85,3.63)
layer2=c(2.03,1.35,2.92,0.31,1.54,0.38,0.54,1.02,1.43,0.96,1.59,1.51,6.68,3.4,8.51,0.66,0.45,0.78,0.07,0.5,0.71,0.7,0.5,1.04,0.09,0,0.86,0.29,1.97,0.4,1.76,2.81,2.47,0,0.62,0.11)
layer3=c(1.29,1.12,1.81,0,0.63,0.25,0.05,0.53,0.38,0,0.36,0,0.23,0.47,0.14,0,0.32,0.02,0,0,0,0,0,0.07,0.67,1.32,0.94,0.18,0.68,0.15,0,1.65,0,1.76,0.34,2.93,0.12,1.26,0.26)
boxplot(layer1,layer2,layer3)

labels = c( replicate(length( layer3 ),"layer3"), replicate(length( layer2 ),"layer2"), replicate(length( layer1 ),"layer1") )
df = data.frame( values=c(layer3,layer2,layer1), labels)


boxplot(values ~ labels, data = df, outpch = NA, ylim=c(0,4)) 
stripchart(values ~ labels, data = df, 
            vertical = TRUE, method = "jitter", 
            pch = 21, col = "maroon", bg = "bisque", 
            add = TRUE) 


> mean(layer1)
[1] 1.102857
> mean(layer2)
[1] 1.415556
> mean(layer3)
[1] 0.5110256
