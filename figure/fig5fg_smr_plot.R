source("input/fig5fg_smr_plot/plot_SMR.r")

if (! dir.exists("output/fig5fg_smr_plot/")) dir.create("output/fig5fg_smr_plot/", recursive = T)

# Fig. 5f: Brain SMR plot
SMR_brain = ReadSMRData("input/fig5fg_smr_plot/Brain_SMR.txt")

pdf("output/fig5fg_smr_plot/fig5f.SMR_brain.pdf")
SMREffectPlot(data=SMR_brain, trait_name="AD")
dev.off()

# Fig. 5g Neuron SMR plot
# note: For optimal visualization in the figure, the axis scale of this plot has been modified from the original PDF output.
SMR_neuron = ReadSMRData("input/fig5fg_smr_plot/Neuron_SMR.txt")

pdf("output/fig5fg_smr_plot/fig5g.SMR_neuron.pdf")
SMREffectPlot(data=SMR_neuron, trait_name="AD") 
dev.off()