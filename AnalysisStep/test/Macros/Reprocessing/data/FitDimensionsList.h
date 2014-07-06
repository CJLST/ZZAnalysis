/*
 * kFitTotalD = 2 x TotalD for each g_i... Thus, kFitTotalD has to be even.
 * Implementation uses kFitTotalD/2 and produces 2 templates for each g_i separately.
 * Template creation requires upon histogram initialization that the second dimensions (D2) to be the interference terms.
 */

const int kFitTotalD = 6;
char* strFitDim[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"D_g1_vs_g2_phi0",
	"D_g1_vs_g4_phi0",
	"D_g1_vs_g2_phi0",
	"D_g2intPdf_phi0",
	"D_g4intPdf_phi0",
	"D_g1_vs_g4_phi0"
};
char* str1DFit_title[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"T_1D_D_g1_vs_g2_phi0",
	"T_1D_D_g1_vs_g4_phi0",
	"T_1D_D_g1_vs_g2_phi0_copy",
	"T_1D_D_g2intPdf_phi0",
	"T_1D_D_g4intPdf_phi0",
	"T_1D_D_g1_vs_g4_phi0_copy"
};
char* str2DFit_title[kFitTotalD/2] = { // Dimension should not exceed kFitTotalD/2!!!
	"T_2D_g1_vs_g2_phi0",
	"T_2D_g1_vs_g4_phi0",
	"T_2D_Dg2_vs_Dg4"
};
char* str3DFit_title[2] = { // Dimension should not exceed 2!!!
	"T_3D_Dg2Dg4Dg2int",
	"T_3D_Dg2Dg4Dg4int"
};

char* strFitDim_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"D_{0^{+}_{h}}",
	"D_{0^{-}}",
	"D_{0^{+}_{h}}",
	"D_{Int}",
	"D_{CP}",
	"D_{CP}"
};
char* str1DFit_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"1D Fit to D_{0^{+}_{h}}",
	"1D Fit to D_{0^{-}}",
	"1D Fit to D_{0^{+}_{h}}",
	"1D Fit to D_{Int}",
	"1D Fit to D_{CP}",
	"1D Fit to D_{CP}"
};
char* str2DFit_label[kFitTotalD/2] = { // Dimension should not exceed kFitTotalD/2!!!
	"2D Fit to D_{0^{+}_{h}} and D_{Int}",
	"2D Fit to D_{0^{-}} and D_{CP}",
	"2D Fit to D_{0^{+}_{h}} and D_{0^{-}}"
};
char* str3DFit_label[2] = { // Dimension should not exceed 2!!!
	"3D Fit to D_{0^{+}_{h}}, D_{0^{-}} and D_{Int}",
	"3D Fit to D_{0^{+}_{h}}, D_{0^{-}} and D_{CP}"
};


const int kFitTotalExtraD = 6;
char* strFitExtraDim[kFitTotalExtraD] = { // Dimension should not exceed kFitTotalExtraD!!!
	"D_g1_vs_g2_phi0",
	"D_g1_vs_g4_phi0",
/* Start extra D templates from here, at strFitExtraDim[t=3], everything above is redundant. */
	"D_g1Q2_phi0",
	"D_g2intPdf_phi90",
	"D_g4intPdf_phi90",
	"D_g1Q2intPdf_phi0"
};
char* str1DExtraFit_title[kFitTotalExtraD] = { // Dimension should not exceed kFitTotalExtraD!!!
	"T_1D_D_g1_vs_g2_phi0",
	"T_1D_D_g1_vs_g4_phi0",
	"T_1D_D_g1Q2_phi0",
	"T_1D_D_g2intPdf_phi90",
	"T_1D_D_g4intPdf_phi90",
	"T_1D_D_g1Q2intPdf_phi0"
};
char* str2DExtraFit_title[kFitTotalExtraD/2] = { // Dimension should not exceed kFitTotalExtraD/2!!!
	"T_2D_g1_vs_g2_withperp",
	"T_2D_g1_vs_g4_withperp",
	"T_2D_g1_vs_g1Q2"
};
char* str3DExtraFit_title[2] = { // Dimension should not exceed 2!!!
	"T_3D_Dg2Dg4Dg2intPerp",
	"T_3D_Dg2Dg4Dg4intPerp"
};

char* strFitExtraDim_label[kFitTotalExtraD] = { // Dimension should not exceed kFitTotalextraD!!!
	"D_{0^{+}_{h}}",
	"D_{0^{-}}",
	"D_{#Lambda1}",
	"D_{Int_{#perp}}",
	"D_{CP_{#perp}}",
	"D_{#Lambda1 Int}"
};
char* str1DExtraFit_label[kFitTotalExtraD] = { // Dimension should not exceed kFitTotalExtraD!!!
	"1D Fit to D_{0^{+}_{h}}",
	"1D Fit to D_{0^{-}}",
	"1D Fit to D_{#Lambda1}",
	"1D Fit to D_{Int_{#perp}}",
	"1D Fit to D_{CP_{#perp}}",
	"1D Fit to D_{#Lambda1 Int}"
};
char* str2DExtraFit_label[kFitTotalExtraD/2] = { // Dimension should not exceed kFitTotalExtraD/2!!!
	"2D Fit to D_{0^{+}_{h}} and D_{Int_{#perp}}",
	"2D Fit to D_{0^{-}} and D_{CP_{#perp}}",
	"2D Fit to D_{#Lambda1} and D_{#Lambda1 Int _{#perp}}"
};
char* str3DExtraFit_label[2] = { // Dimension should not exceed 2!!!
	"3D Fit to D_{0^{+}_{h}}, D_{0^{-}} and D_{Int_{#perp}}",
	"3D Fit to D_{0^{+}_{h}}, D_{0^{-}} and D_{CP_{#perp}}"
};

char strDBkg[]="D_bkg";
char* str2DBkgFit_title[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"T_2DBkg_D_g1_vs_g2_phi0",
	"T_2DBkg_D_g1_vs_g4_phi0",
	"T_2DBkg_D_g1Q2_phi0",
	"T_2DBkg_D_g2intPdf_phi0",
	"T_2DBkg_D_g4intPdf_phi0",
	"T_2DBkg_D_g1Q2intPdf_phi0"
};
char* str2DBkgFit_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"2D Fit to D_{0^{+}_{h}} and D_{Bkg}",
	"2D Fit to D_{0^{-}} and D_{Bkg}",
	"2D Fit to D_{#Lambda1} and D_{Bkg}",
	"2D Fit to D_{Int} and D_{Bkg}",
	"2D Fit to D_{CP} and D_{Bkg}",
	"2D Fit to D_{#Lambda1 Int} and D_{Bkg}"
};
