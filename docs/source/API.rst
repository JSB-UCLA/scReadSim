API
===

Import scReadSim as::

   import scReadSim


Utility
~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.Utility.CallPeak
   scReadSim.Utility.ExtractBAMCoverage
   scReadSim.Utility.scATAC_CreateFeatureSets
   scReadSim.Utility.scRNA_CreateFeatureSets
   scReadSim.Utility.countmat_mainloop
   scReadSim.Utility.scATAC_bam2countmat_paral
   scReadSim.Utility.scRNA_UMIcountmat_mainloop
   scReadSim.Utility.scRNA_bam2countmat_paral
   scReadSim.Utility.bam2countmat_INPUT
   scReadSim.Utility.find_nearest
   scReadSim.Utility.match_peak
   scReadSim.Utility.ComplementFeature


scATAC_GenerateBAM
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.scATAC_GenerateBAM.flatten
   scReadSim.scATAC_GenerateBAM.cellbarcode_generator
   scReadSim.scATAC_GenerateBAM.scATAC_INPUT_PerTruePeakEdition
   scReadSim.scATAC_GenerateBAM.scATAC_PerTruePeakEdition
   scReadSim.scATAC_GenerateBAM.generateBAMcoord_mainloop
   scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord_paral
   scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT
   scReadSim.scATAC_GenerateBAM.scATAC_CombineBED
   scReadSim.scATAC_GenerateBAM.scATAC_BED2FASTQ
   scReadSim.scATAC_GenerateBAM.AlignSyntheticBam_Pair
   scReadSim.scATAC_GenerateBAM.ErrorBase
   scReadSim.scATAC_GenerateBAM.ErroneousRead
   scReadSim.scATAC_GenerateBAM.SubstiError_Pair
   scReadSim.scATAC_GenerateBAM.scATAC_ErrorBase


scRNA_GenerateBAM
~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.scRNA_GenerateBAM.flatten
   scReadSim.scRNA_GenerateBAM.cellbarcode_generator
   scReadSim.scRNA_GenerateBAM.scRNA_SampleSyntheticReads
   scReadSim.scRNA_GenerateBAM.scRNA_PerTruePeakEdition_read
   scReadSim.scRNA_GenerateBAM.scRNA_PerTruePeakEdition_UMI
   scReadSim.scRNA_GenerateBAM.generateBAMcoord_read_mainloop
   scReadSim.scRNA_GenerateBAM.generateBAMcoord_UMI_mainloop
   scReadSim.scRNA_GenerateBAM.scRNA_GenerateBAMCoord_paral
   scReadSim.scRNA_GenerateBAM.scRNA_CombineBED
   scReadSim.scRNA_GenerateBAM.scRNA_BED2FASTQ
   scReadSim.scRNA_GenerateBAM.AlignSyntheticBam_Single
   scReadSim.scRNA_GenerateBAM.ErrorBase
   scReadSim.scRNA_GenerateBAM.ErroneousRead
   scReadSim.scRNA_GenerateBAM.SubstiError
   scReadSim.scRNA_GenerateBAM.scRNA_ErrorBase


GenerateSyntheticCount
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.GenerateSyntheticCount.scATAC_GenerateSyntheticCount
   scReadSim.GenerateSyntheticCount.scRNA_GenerateSyntheticCount


   


