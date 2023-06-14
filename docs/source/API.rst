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
   scReadSim.Utility.scATAC_bam2countmat_OutputPeak
   scReadSim.Utility.find_nearest_peak
   scReadSim.Utility.find_nearest_nonpeak
   scReadSim.Utility.match_peak
   scReadSim.Utility.match_nonpeak
   scReadSim.Utility.bam2MarginalCount
   scReadSim.Utility.FeatureMapping


scATAC_GenerateBAM
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.scATAC_GenerateBAM.flatten
   scReadSim.scATAC_GenerateBAM.cellbarcode_generator
   scReadSim.scATAC_GenerateBAM.find_leftnearest_nonpeak
   scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord
   scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak
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
   scReadSim.scRNA_GenerateBAM.scRNA_GenerateBAMCoord
   scReadSim.scRNA_GenerateBAM.scRNA_CombineBED
   scReadSim.scRNA_GenerateBAM.scRNA_BED2FASTQ
   scReadSim.scRNA_GenerateBAM.AlignSyntheticBam_Single
   scReadSim.scRNA_GenerateBAM.ErrorBase
   scReadSim.scRNA_GenerateBAM.ErroneousRead
   scReadSim.scRNA_GenerateBAM.SubstiError
   scReadSim.scRNA_GenerateBAM.scRNA_ErrorBase

DoubletDetection
~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.DoubletDetection.detectDoublet


GenerateSyntheticCount
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.GenerateSyntheticCount.scATAC_GenerateSyntheticCount
   scReadSim.GenerateSyntheticCount.scRNA_GenerateSyntheticCount


GenerateSyntheticCount_MultiOmics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.GenerateSyntheticCount_MultiOmics.scMultiOmics_GenerateSyntheticCount

   


