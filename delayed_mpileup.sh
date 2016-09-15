SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
while kill -0 36485; do
    sleep 60
done
nohup bash $SCRIPTS/create_mpileup.sh -l fem_vie2011_pileup.log -p fem_vie2011 pI25=../PopI/ColdShock25/BGI_34_indx2_CGATGTAT/BGI_34_indx2_CGATGTAT.picnodupl.filtered.mq20_chrfilt.bam pII25=../PopII/ColdShock25/BGI_34_indx4_TGACCAAT/BGI_34_indx4.picnodupl.filtered.mq20_chrfilt.bam pIII25=./PopIII/ColdShock25/BGI_34_indx5_ACAGTGAT/BGI_34_indx5_ACAGTGAT.picnodupl.filtered.mq20_chrfilt.bam pIp=../PopI/Pupation/BGI_2_hl_merged.bam pIIp=../PopII/Pupation/BGI_21a_3a_merged.bam pIIIp=../PopIII/Pupation/BGI_22a_23a_merged.bam &
exit 0