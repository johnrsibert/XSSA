#/bin/bash
# draw map of US EEZ around Hawaii showing Main Hawaiian Islands and
# MFCL regions 2 and 4

MAP=HI_regions
echo $MAP
EPS=$MAP.eps
echo $EPS
rm -v $EPS
rm -v $MAP.png

REGION=-R175/210/15/35
echo $REGION
PROJECTION=-JM6i
echo $PROJECTION

rm -fv .gmt*
gmtset ANOT_FONT_SIZE_PRIMARY 12 BASEMAP_TYPE fancy PAGE_ORIENTATION portrait HEADER_FONT_SIZE 16 HEADER_FONT 1 MEASURE_UNIT inch PLOT_DEGREE_FORMAT dddF

# shift from west longitude to east longitude
rm -fv junk*.txt
awk '{ printf("%lf %lf\n", 360-$1, $2)}' hi_eez.txt > junk.txt 
# dig out MHI outline
awk '{ if ($1 > 198.67) printf("%lf %lf\n", $1, $2)}' junk.txt > junk2.txt

# draw basemap and border
psbasemap -B10f5 $REGION $PROJECTION  -V    -K > $EPS 
# draw shaded MHI
psxy $REGION $PROJECTION junk2.txt  -Gp300/7  -V -O -K >> $EPS
# draw US EEZ
psxy $REGION $PROJECTION junk.txt  -W2p  -V -O -K >> $EPS

# draw MHI boundary
psxy $REGION $PROJECTION -A -W1p -V -O -K  << EOF >> $EPS
198.67 26.31
198.67 18.33
EOF

# draw MFCL regions 2 & 4 boundary
MFCL=black
psxy $REGION $PROJECTION -A -W1p,$MFCL -V -O -K  << EOF >> $EPS
210 20
175 20
EOF


# label MFCL regions 2 & 4
echo '185.00  21.50 10 0 0 6 MFCL 2'| pstext -G$MFCL  $REGION $PROJECTION -V -O -K >> $EPS
echo '185.00  18.50 10 0 0 6 MFCL 4'| pstext -G$MFCL  $REGION $PROJECTION -V -O -K >> $EPS

#plot high resolution coast map
pscoast $REGION $PROJECTION -Dh  -G237/214/151  -W1,154/139/98 -V -O    >> $EPS

#gv $EPS &
ps2raster -V -A -Tg -E600 $EPS
eog $MAP.png &
