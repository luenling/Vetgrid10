get_FDR_cutoff = function(obs_pvals,null_pvals,ranks=1:100,fdr=0.05) { 
	obs_pvals=sort(obs_pvals)
	n_ratio=length(obs_pvals)/length(null_pvals)
	result= data.frame(rank = numeric(0), obs_P = numeric(0), FDR = numeric(0))
	for (r in ranks)
		{
		null_pvals_smaller = n_ratio* sum(null_pvals < obs_pvals[r])
		result[r,]=c(r,obs_pvals[r], null_pvals_smaller/r)
		if (result[r,3] >= fdr && r) {break}
		}
	return(result)
	}
# comparing shuffled to observed p-Values
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/CMH")
cmh_joined=read.table("fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.gwas",header=T)

pV_shuf = scan("shuffles/selected_shuffle_ps.dat")
quartz()
çabline(a=0,b=1,col="red")
dev.print(png,width=800,file="qq_plot_vie_ita_joined_shuffled_14shuff.png")
obs_pvals=cmh_joined$P
obs_pvals=sort(obs_pvals)
ranks = 200:250
n_ratio=length(obs_pvals)/length(pV_shuf)
for (r in ranks)
	{
	pV_shuf_smaller = n_ratio* sum(pV_shuf < obs_pvals[r])
	cat ("r=", r, obs_pvals[r], pV_shuf_smaller/r, "\n")
	}

setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/")
cmh_ita=read.table("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.gwas",header=T)
pV_shuf=scan("shuffles/selected_shuffle_ps.dat")
nlpV=-1*log10(pV_shuf)
nl_ita=-1*log10(cmh_ita$P)
l_ita=length(nl_ita)
quants = c(seq(1,200), seq(201,1000,100),seq(1001,10000,1000),seq(10001,100000,10000),seq(100001,l_ita,100000))/l_ita
q_ita=quantile(cmh_ita$P,  probs = quants)
q_pV=quantile(pV_shuf,  probs = quants)
plot(-log10(q_pV),-log10(q_ita),xlab="null P (8 shuffles)", ylab="observed P", main="Italy 2011, shuffles with <= 2 pairs")
abline(a=0,b=1,col="red")
dev.print(png,width=800,file="qq_plot_ita_shuffled_8shuff_sparse.png")
qqplot(-1*log10(pV_shuf),-1 * log10(cmh_ita$P),xlab="null P (8 shuffles)", ylab="observed P", main="Italy 2011, shuffles with <= 2 pairs")
abline(a=0,b=1,col="red")
dev.print(png,width=800,file="qq_plot_vie_ita_joined_shuffled_8shuff.png")
results=
