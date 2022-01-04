#for plotting 
#load("D:/PhD_NL/Codes/mutualism/updated_state.Rdata") 
#load("D:/PhD_NL/Codes/mutualism/rates.Rdata")

Mt <- updated_state$Mt
p_status <- updated_state$p_status
a_status <- updated_state$a_status
island_spec_plant <- updated_state$island_spec_plant
island_spec_animal <- updated_state$island_spec_animal
maxplantID <- updated_state$maxplantID
maxanimalID <- updated_state$maxanimalID
stt_table <- updated_state$stt_table

#stt_table: TIme nIp nAp nCp nIa nAa nCa 

par(mfrow=c(2,3)) 
#plot(x,y1,type='o'/'p'/'l', xlim=[], ylim=[], xlab='',lty=1, lwd=2, main="title",axes = TRUE,
#     pch='', cex=, col=red) 
plot(stt_table[,2]~stt_table[,1], xlab='Time', ylab='nIp', main = "nIp_Time")
plot(stt_table[,3]~stt_table[,1], xlab='Time', ylab='nAp',  main = "nAp_Time")
plot(stt_table[,4]~stt_table[,1], xlab='Time', ylab='nCp',  main = "nCp_Time")
plot(stt_table[,5]~stt_table[,1], xlab='Time', ylab='nIa',  main = "nIa_Time")
plot(stt_table[,6]~stt_table[,1], xlab='Time', ylab='nAa',  main = "nAa_Time")
plot(stt_table[,7]~stt_table[,1], xlab='Time', ylab='nCa',  main = "nCa_Time")

# endemic & non-endemic for plant species and animal species
par(mfrow=c(1,2)) 

plot(stt_table[,2]~stt_table[,1], type='l', col = 'red', xlab='Time', main = "plant species")
lines(stt_table[,3]+stt_table[,4]~stt_table[,1],col='blue',xlab='Time')

plot(stt_table[,5]~stt_table[,1], type='l', col = 'red', xlab='Time', main = "animal species")
lines(stt_table[,6]+stt_table[,7]~stt_table[,1],col='blue',xlab='Time')






