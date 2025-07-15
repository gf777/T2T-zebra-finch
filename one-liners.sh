######## extra ########
# basic dotplot
mashmap -r ref.fa -q qry.fa --pi 95 -s 10000 -t 32 -o mashmap.out
perl generateDotPlot png large mashmap.out

# to better visualize sequence direction in bandage remember to:
- load graph into bandageNG
- Select BandageNG-> Preferences
- Under graph appearence -> "Arrowheads in sigle node style" switch to on.
- Draw graph as usual