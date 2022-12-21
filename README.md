# MSstats-tutorial
MSstats tutorial from the May Institute, 2020

Based on the following video;
https://www.youtube.com/watch?v=30a6sDj8Qlo&list=PL2u38g_AG4MG_puKZBT4HawiVLQBE28gD&index=10

Github for tutorial;
https://meenachoi.github.io/MayInstitute2019RstatsMS/

Data from the following;
https://docs.google.com/document/d/1vnMG0XTqYVngpK5EyIQNXi-XDDh7sq72yQnhrnxw83k/edit

Further talks found here;
https://computationalproteomics2020.khoury.northeastern.edu/

And talks from 2022;
https://computationalproteomics.khoury.northeastern.edu/

Google group for troubleshooting `MSstats`;
https://groups.google.com/g/msstats

## Ideas and thoughts
It might be easiest to do basic QC using the artMS package;
https://www.bioconductor.org/packages/devel/bioc/vignettes/artMS/inst/doc/artMS_vignette.html#1_OVERVIEW

## Current issues
Currently I'm struggling to get the function `PDtoMSstatsFormat` to work. Initially I was missing the columns for "ProteinGroupAccessions" and "PrecursorArea", however I renamed the column"Master.Protein.Accessions" from the original dataset. "PrecursorArea" can be renamed from the "Precursor.Area" or "Precursor.Abundance" column, but note that these are only produced when running a quantified search in Proteome Discoverer.
However, I'm still having other issues, being discussed here:
https://groups.google.com/g/msstats/c/fHa3MqRtMss/m/PCLGjy0SCAAJ