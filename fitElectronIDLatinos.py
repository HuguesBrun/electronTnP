import FWCore.ParameterSet.Config as cms
import re, os
### USAGE:
###    cmsRun fitElecID.py <scenario> [ <id> [ <binning1> ... <binningN> ] ]
###
### scenarios:
###   - data_all (default)  
###   - signal_mc

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-elec Mass", "60", "120", "GeV/c^{2}"),
        pt = cms.vstring("elec p_{T}", "0", "10000", "GeV/c"),
        eta    = cms.vstring("elec #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("elec |#eta|", "0", "2.5", ""),
        absSCeta = cms.vstring("elec  |SC #eta|", "0","2.5",""),
        SCeta = cms.vstring("elec SC #eta", "-2.5","2.5",""),
        nVtx   = cms.vstring("Number of vertices", "0", "999", ""),
        weight = cms.vstring("weight","0","1000",""),
        weight_runA = cms.vstring("weight_runA","0","1000",""),
        weight_runB = cms.vstring("weight_runB","0","1000",""),
        weight_runC = cms.vstring("weight_runC","0","1000",""),
        weight_runD = cms.vstring("weight_runD","0","1000",""),
       ),

    Categories = cms.PSet(
                          passTight   = cms.vstring("passTight", "dummy[pass=1,fail=0]"),
                          passLoose   = cms.vstring("passLoose", "dummy[pass=1,fail=0]"),
                          passFO   = cms.vstring("passFO", "dummy[pass=1,fail=0]"),
                          passBDT   = cms.vstring("passBDT", "dummy[pass=1,fail=0]"),
                          passISO   = cms.vstring("passISO", "dummy[pass=1,fail=0]"),
                          passFO_BDT   = cms.vstring("passFO_BDT", "dummy[pass=1,fail=0]"),
                          passFO_ISO   = cms.vstring("passFO_ISO", "dummy[pass=1,fail=0]"),
                          passFO_BDT_ISO   = cms.vstring("passFO_BDT_ISO", "dummy[pass=1,fail=0]"),
                          passPreselec   = cms.vstring("passPreselec", "dummy[pass=1,fail=0]"),
                          passIP   = cms.vstring("passIP", "dummy[pass=1,fail=0]"),
                          passConvs   = cms.vstring("passConvs", "dummy[pass=1,fail=0]"),
                          passFOnoIso   = cms.vstring("passFOnoIso", "dummy[pass=1,fail=0]"),
                          passBDT_ISO   = cms.vstring("passBDT_ISO", "dummy[pass=1,fail=0]"),
                          passAllNoIsoDet   = cms.vstring("passAllNoIsoDet", "dummy[pass=1,fail=0]"),
                          passNM1IP   = cms.vstring("passNM1IP", "dummy[pass=1,fail=0]"),
                          passNM1convs   = cms.vstring("passNM1convs", "dummy[pass=1,fail=0]"),
                          passNM1presel   = cms.vstring("passNM1presel", "dummy[pass=1,fail=0]"),
                          trigSingle = cms.vstring("trigSingle", "dummy[pass=1,fail=0]"),
                          trigDoubleLeg0 = cms.vstring("trigDoubleLeg0", "dummy[pass=1,fail=0]"),
                          trigDoubleLeg1 = cms.vstring("trigDoubleLeg1", "dummy[pass=1,fail=0]"),
                          eventMatched = cms.vstring("eventMatched", "dummy[pass=1,fail=0]"),
                          isSameSign = cms.vstring("isSameSign","dummy[pass=1,fail=0]")
    ),
    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpoMin70 = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,3,10])",
            "SUM::signal(vFrac[0.8,0.5,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0.7,1]",
            "signalFractionInPassing[0.9]"
        ),
        CBBreitWignerPlusExponentialBackground = cms.vstring(
              "mean[90,85,95]",
              "scale[1,0.,2.]",
              "sigma[3,1,10]",
              "sigmaF[3,1,10]",
              "CBShape::cbs(mass, scale, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
              "CBShape::cbsF(mass, scale, sigmaF, alphaF[3., 0.5, 5.], n[1, 0., 100.])",
              "RooBreitWigner::vs(mass, mean, sigma)",
              "RooBreitWigner::vsF(mass, mean, sigmaF)",
              "FCONV::signalPass(mass,vs,cbs)",
              "FCONV::signalFail(mass,vsF,cbsF)",
              "Exponential::backgroundPass(mass, lp[0,-5,5])",
              "Exponential::backgroundFail(mass, lf[0,-5,5])",
              "efficiency[0.9,0,1]", 
              "signalFractionInPassing[0.9,0.5,1]"
       ),
        CBBreitWignerPlusPolyBackground = cms.vstring(
              "mean[90,85,95]",
              "scale[1,0.05,2.]",
              "sigma[3,1,10]",
              "sigmaF[3,1,10]",
              "CBShape::cbs(mass, scale, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
              "CBShape::cbsF(mass, scale, sigmaF, alphaF[3., 0.5, 5.], n[1, 0., 100.])",
              "RooBreitWigner::vs(mass, mean, sigma)",
              "RooBreitWigner::vsF(mass, mean, sigmaF)",
              "FCONV::signalPass(mass,vs,cbs)",
              "FCONV::signalFail(mass,vsF,cbsF)",
              "RooLandau::backgroundPass(mass, Lmp[100,90,105],wp[1,0,10])",
              "RooLandau::backgroundFail(mass, Lmf[100,90,105],wf[1,0,10])",
              "efficiency[0.9,0,1]",
              "signalFractionInPassing[0.9,0.5,1]"
       ),
       SignalFromMC_Landau = cms.vstring(
              "mean[90,85,95]",
              "scale[1,0.05,2.]",
              "sigma[3,1,10]",
              "sigmaF[3,1,10]",
              "CBShape::cbs(mass, scale, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
              "CBShape::cbsF(mass, scale, sigmaF, alphaF[3., 0.5, 5.], n[1, 0., 100.])",
              "RooBreitWigner::vs(mass, mean, sigma)",
              "RooBreitWigner::vsF(mass, mean, sigmaF)",
              "FCONV::convPass(mass,vs,cbs)",
              "FCONV::convFail(mass,vsF,cbsF)",
              "RooLandau::landp(mass, Lmp[100,90,105],wp[1,0,10])",
              "RooLandau::landf(mass, Lmf[100,90,105],wf[1,0,10])",
	          "SUM::signalPass(vFrac[0.8,0.5,1]*convPass, landp)",	
              "SUM::signalFail(vFrac[0.8,0.5,1]*convFail, landf)",	
              "Exponential::backgroundPass(mass, lp[0])",
              "Exponential::backgroundFail(mass, lf[0])",
              "efficiency[0.9,0,1]",
              "signalFractionInPassing[1]"
       ), 
       SignalFromMC_BadZone = cms.vstring(
              "mean[90,85,95]",
              "Voigtian::signal1p(mass, mean1p[90,80,100], widthp[2.495], sigma1p[2,1,3])",
              "Voigtian::signal2p(mass, mean2p[90,80,100], widthp,        sigma2p[4,2,10])",
              "SUM::Sumsignalp(vPropP[0.8,0,1]*signal1p, signal2p)",
              "Voigtian::signal1f(mass, mean1f[90,80,100], widthf[2.495], sigma1f[2,1,3])",
              "Voigtian::signal2f(mass, mean2f[90,80,100], widthf,        sigma2f[4,2,10])",
              "SUM::Sumsignalf(vPropF[0.8,0,1]*signal1f, signal2f)",
              "Exponential::expP(mass, lep[0,-5,5])",
              "Exponential::expF(mass, lef[0,-5,5])",
              "SUM::signalPass(vFrac[0.8,0.5,1]*Sumsignalp, expP)",	
              "SUM::signalFail(vFrac[0.8,0.5,1]*Sumsignalf, expF)",
              "Exponential::backgroundPass(mass, lp[0])",
              "Exponential::backgroundFail(mass, lf[0])",
              "efficiency[0.9,0,1]",
              "signalFractionInPassing[1]"
        ), 
        SignalFromMC = cms.vstring(
              "mean[90,85,95]",
              "scale[1,0.,2.]",
              "sigma[3,1,10]",
              "sigmaF[3,1,10]",
              "CBShape::cbs(mass, scale, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
              "CBShape::cbsF(mass, scale, sigmaF, alphaF[3., 0.5, 5.], n[1, 0., 100.])",
              "RooBreitWigner::vs(mass, mean, sigma)",
              "RooBreitWigner::vsF(mass, mean, sigmaF)",
              "FCONV::convPass(mass,vs,cbs)",
              "FCONV::convFail(mass,vsF,cbsF)",
              "Exponential::expP(mass, lep[0,-5,5])",
              "Exponential::expF(mass, lef[0,-5,5])",
	          "SUM::signalPass(vFrac[0.8,0.5,1]*convPass, expP)",	
              "SUM::signalFail(vFrac[0.8,0.5,1]*convFail, expF)",	
              "Exponential::backgroundPass(mass, lp[0])",
              "Exponential::backgroundFail(mass, lf[0])",
              "efficiency[0.9,0,1]",
              "signalFractionInPassing[1]"
       ),
                   SignalFromMCMyBgShape = cms.vstring(
                                                      "mean[90,85,95]",
                                                      "sigma[3,1,10]",
                                           #  "meanF[90,85,95]",
                                           #          "sigmaF[3,1,10]",
                                           #            "Exponential::signal(mass, lp[1,1,5])",
                                           "CBShape::signal(mass, mean, sigma, alpha[3., 0.5, 5.], n[1, 0., 100.])",
                                           #      "CBShape::signal(mass, meanF, sigmaF, alphaF[3., 0.5, 5.], n[1, 0., 100.])",
                                           #          "Exponential::expPart(mass, lf[-1,-5,0])",
                                           "RooBernstein::backgroundFail(mass, {a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50],a4[1,0,50]})",
                                           "RooBernstein::backgroundPass(mass, {a0p[10,0,50],a1p[1,0,50],a2p[1,0,50],a3p[1,0,50],a4p[1,0,50]})",
                                           # "SUM::backgroundFail(vFrac[0.8,0.5,1]*theCBPart,bernPart)",
                                                      "efficiency[0.9,0,1]",
                                                      "signalFractionInPassing[0.2,0,1]"
                                                      ), 
            ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),

    Efficiencies = cms.PSet(), # will be filled later
)





TRIGGER = cms.PSet(tag_Mu24 = cms.vstring("pass"))
if "mc" in scenario or "39X" in scenario or "38X" in scenario:
    TRIGGER = cms.PSet(tag_Mu15 = cms.vstring("pass"), tag_pt = cms.vdouble(24.,9999.))

PT_ETA_BINS = cms.PSet(
    pt     = cms.vdouble(  15, 20, 25, 50, 150 ),
    absSCeta = cms.vdouble(  0.0, 1.4442, 1.556, 2.5),
)


PT_ETA_BINS_XCHECK = cms.PSet(
                              pt  = cms.vdouble(20, 30, 40, 60, 100),
                              abseta = cms.vdouble(0, 1.479, 2.5),
                              )


VTX_BINS  = cms.PSet(
    pt     = cms.vdouble(  10, 30, 50 ),
    abseta = cms.vdouble(  0.0, 2.4),
	nVtx = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5, 32.5, 34.5),

)

PT_ETA_BINS_POG = cms.PSet(
                                  pt = cms.vdouble(10,15,20,30,40,50,200),
                                  absSCeta = cms.vdouble(0, 0.8, 1.4442, 1.556, 2, 2.5),
                                  isSameSign = cms.vstring("fail"),
                                  )


if "mc" in scenario:
    setattr(PT_ETA_BINS_POG, "eventMatched",cms.vstring("pass"))

if "NM1ISO" in scenario:
    setattr(PT_ETA_BINS_POG, "passFO_BDT",cms.vstring("pass"))
if "NM1ID" in scenario:
    setattr(PT_ETA_BINS_POG, "passFO_ISO",cms.vstring("pass"))

PREFIX="file:/afs/cern.ch/work/h/hbrun/TnPmyProdTer/"
#PREFIX="file:/tmp/hbrun/"
process.TnP_ElecID = Template.clone(
    InputFileNames = cms.vstring("TnPTree.root"),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("TnP_ElecID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)
if "mc" in scenario:
    if "RunA" in scenario: process.TnP_ElecID.WeightVariable = cms.string("weight_runA")
    elif "RunB" in scenario: process.TnP_ElecID.WeightVariable = cms.string("weight_runB")
    elif "RunC" in scenario: process.TnP_ElecID.WeightVariable = cms.string("weight_runC")
    elif "RunD" in scenario: process.TnP_ElecID.WeightVariable = cms.string("weight_runD")

if "data" in scenario:
    if   "v1" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpZ_HWWid2012.root" ]
    elif "v2" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"tnpZ_2011A_v2_GOLDEN.root" ]
    elif "huguesTestRunA" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"TnP_GoodShapeelectrons_runA.root"]
    elif "huguesTestRunB" in scenario: process.TnP_ElecID.InputFileNames = [PREFIX+"TnP_GoodShapeelectrons_runB_part0.root",PREFIX+"TnP_GoodShapeelectrons_runB_part1.root",PREFIX+"TnP_GoodShapeelectrons_runB_part2.root"]
    elif "huguesTestRunC" in scenario: process.TnP_ElecID.InputFileNames = [PREFIX+"TnP_GoodShapeelectrons_runCv1_RERECO.root",PREFIX+"TnP_GoodShapeelectrons_runCv2_part0.root",PREFIX+"TnP_GoodShapeelectrons_runCv2_part1.root",PREFIX+"TnP_GoodShapeelectrons_runCv2_part2.root",PREFIX+"TnP_GoodShapeelectrons_runCv2_part3.root",PREFIX+"TnP_GoodShapeelectrons_runCv2_part4.root"]
    elif "huguesTestRunD" in scenario: process.TnP_ElecID.InputFileNames = [ PREFIX+"TnP_GoodShapeelectrons_runD_part1.root",PREFIX+"TnP_GoodShapeelectrons_runD_part2.root",PREFIX+"TnP_GoodShapeelectrons_runD_part3.root"]
    else: print "coucou c encore moi"
if "mc" in scenario:
    process.TnP_ElecID.InputFileNames = [PREFIX+"TnP_GoodShapeDYJet_ter3.root" ]

if "tag35" in scenario:
    process.TnP_ElecID.Variables.tag_pt[1]='35'

print "les fichiers que l'on va utiliser = ", process.TnP_ElecID.InputFileNames

#IDS = ["TOGCPFTIPMVA_from_TrackerOrGlobal"]
#IDS = [ "passPogLoose","passPogTight","passTrig","passTrigPlusID","passIdTrigIso","passIsoTrig","passID2012_from_passTrig","passISO2012_from_passTrig","passISO2012_from_passID2012","passID2012_from_passISO2012"]
#IDS = ["dPhiCut","dEtaCut","IPcut","isoECALCut","isoHCAL03Cut","isoTrackerCut","HoEcut","sigmaIEtaIetaCut","passConvR"]
#test with no ECAL cutIDS = ["passFOselection_noEcal","passFOselection","passAll_noEcal_from_passISO2012","passElec_FO_ISO_ID_from_passISO2012","passAll_noEcal_from_passISO2012","passElec_FO_ISO_ID_from_passISO2012"]
#IDS = ["passTight","passFO_BDT_ISO"]
IDS = ["passFO_BDT_ISO"]#,"passTight"]
if "_FO" in scenario: IDS = ["passFO"]


#,"passElec_FO","passElec_FO_ID","passElec_FO_ISO_ID_from_passISO2012","passElec_FO_ISO_ID_from_passElec_FO_ID","passElec_FO_ISO_ID_from_passElec_FO","passPogTight","passPogLoose"]
#ALLBINS = [("ptXcheck",PT_ETA_BINS_XCHECK)]


TemplateSignal_Landau = cms.vstring(
                                    "scaleTp[1,0.9,1.1]",
                                    "scaleTf[1,0.9,1.1]",
                                    "largerResPass[1,0.,2.]",
                                    "largerResFail[1,0.,2.]",
                                    "expr::NewMean1p('mean*scaleTp',mean,scaleTp)",
                                    "expr::NewMean1f('mean*scaleTf',mean,scaleTf)",
                                    "expr::NewSigma1p('sigma*largerResPass',sigma,largerResPass)",
                                    "expr::NewSigma1f('sigmaF*largerResFail',sigmaF,largerResFail)",
                                    "CBShape::cbs(mass, scale, NewSigma1p, alpha, n)",
                                    "CBShape::cbsF(mass, scale, NewSigma1f, alphaF, n)",
                                    "RooBreitWigner::vs(mass, NewMean1p, NewSigma1p)",
                                    "RooBreitWigner::vsF(mass, NewMean1f, NewSigma1f)",
                                    "FCONV::convPass(mass,vs,cbs)",
                                    "FCONV::convFail(mass,vsF,cbsF)",
                                    "RooLandau::landp(mass, Lmp,wp)",
                                    "RooLandau::landf(mass, Lmf,wf)",
                                    "SUM::signalPass(vFrac*convPass, landp)",
                                    "SUM::signalFail(vFrac*convFail, landf)",
                                    "RooBernstein::backgroundPass(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                                    "RooBernstein::backgroundFail(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                                    # "Exponential::backgroundPass(mass, lp[0,-5,5])",
                                    #"Exponential::backgroundFail(mass, lf[0,-5,5])",
                                    "efficiency[0.9,0,1]",
                                    "signalFractionInPassing[1]"
                                      )
TemplateSignal = cms.vstring(
                             "scaleTp[1,0.9,1.1]",
                             "scaleTf[1,0.9,1.1]",
                             "largerResPass[1,0.,2.]",
                             "largerResFail[1,0.,2.]",
                             "expr::NewMean1p('mean*scaleTp',mean,scaleTp)",
                             "expr::NewMean1f('mean*scaleTf',mean,scaleTf)",
                             "expr::NewSigma1p('sigma*largerResPass',sigma,largerResPass)",
                             "expr::NewSigma1f('sigmaF*largerResFail',sigmaF,largerResFail)",
                             "CBShape::cbs(mass, scale, NewSigma1p, alpha, n)", #[3., 0.5, 5.]
                             "CBShape::cbsF(mass, scale, NewSigma1f, alphaF, n)",
                             "RooBreitWigner::vs(mass, NewMean1p, NewSigma1p)",
                             "RooBreitWigner::vsF(mass, NewMean1f, NewSigma1f)",
                             "FCONV::convPass(mass,vs,cbs)",
                             "FCONV::convFail(mass,vsF,cbsF)",
                             "Exponential::expP(mass, lep)",
                             "Exponential::expF(mass, lef)",
                             "SUM::signalPass(vFrac*convPass, expP)",
                             "SUM::signalFail(vFrac*convFail, expF)",
                             #"Exponential::backgroundFail(mass, lf[0,-5,5])",
                             #"Exponential::backgroundPass(mass, lp[0,-8,8])",
                             "RooBernstein::backgroundPass(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                             "RooBernstein::backgroundFail(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                             "efficiency[0.9,0,1]",
                             "signalFractionInPassing[0.9,0.5,1]"
                             ) 




TemplateSignal_BadZone = cms.vstring(
                                     "scaleTp[1,0.9,1.1]",
                                     "scaleTf[1,0.9,1.1]",
                                     "largerResPass[1,0.,2.]",
                                       "largerResFail[1,0.,2.]",
                                       "expr::NewMean1p('mean1p*scaleTp',mean1p,scaleTp)",
                                       "expr::NewMean2p('mean2p*scaleTp',mean2p,scaleTp)",
                                       "expr::NewMean1f('mean1f*scaleTf',mean1f,scaleTf)",
                                       "expr::NewMean2f('mean2f*scaleTf',mean2f,scaleTf)",
                                       "expr::NewSigma1p('sigma1p*largerResPass',sigma1p,largerResPass)",
                                       "expr::NewSigma2p('sigma2p*largerResPass',sigma2p,largerResPass)",
                                       "expr::NewSigma1f('sigma1f*largerResFail',sigma1f,largerResFail)",
                                       "expr::NewSigma2f('sigma2f*largerResFail',sigma2f,largerResFail)",
                                       "Voigtian::signal1p(mass, NewMean1p, widthp[2.495], NewSigma1p)",
                                       "Voigtian::signal2p(mass, NewMean2p, widthp,        NewSigma2p)",
                                       "SUM::Sumsignalp(vPropP*signal1p, signal2p)",
                                       "Voigtian::signal1f(mass, NewMean1f, widthf[2.495], NewSigma1f)",
                                       "Voigtian::signal2f(mass, NewMean2f, widthf,        NewSigma2f)",
                                       "SUM::Sumsignalf(vPropF*signal1f, signal2f)",
                                       "Exponential::expP(mass, lep)",
                                       "Exponential::expF(mass, lef)",
                                       "SUM::signalPass(vFrac*Sumsignalp, expP)",	
                                       "SUM::signalFail(vFrac*Sumsignalf, expF)",
                                     #"Exponential::backgroundPass(mass, lp[-7,-8,8])",
                                     #  "Exponential::backgroundFail(mass, lf[-7,-8,8])",
                                       "RooBernstein::backgroundPass(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                                       "RooBernstein::backgroundFail(mass,{a0[10,0,50],a1[1,0,50],a2[1,0,50],a3[1,0,50]})",
                                       "efficiency[0.9,0,1]",
                                       "signalFractionInPassing[0.9,0.5,1]"
                                       )
###### now built the pdfs !!!! 
if "_FO" in scenario: theCatPart="FO"
if "_NM1ISO" in scenario: theCatPart="NM1ISO"
if "_NM1ID" in scenario: theCatPart="NM1ID"


if "data" in scenario and not "forBg" in scenario:
    for ID in IDS:
        file = open(ID+"_pt"+theCatPart+".info","r")
        files = file.readlines()
        file.close()
        for line in files:
            ligneSplitted = re.split(" ",line)
            if (("pt_bin5" in ligneSplitted[1]) and not ("absSCeta_bin2" in ligneSplitted[0])):
                localPDF = cms.vstring()
                localPDF.extend(["Lmf["+str(ligneSplitted[2])+"]"])
                localPDF.extend(["Lmp["+str(ligneSplitted[3])+"]"])
                localPDF.extend(["alpha["+str(ligneSplitted[4])+"]"])
                localPDF.extend(["alphaF["+str(ligneSplitted[5])+"]"])
                localPDF.extend(["mean["+str(ligneSplitted[9])+"]"])
                localPDF.extend(["n["+str(ligneSplitted[10])+"]"])
                localPDF.extend(["scale["+str(ligneSplitted[12])+"]"])
                localPDF.extend(["sigma["+str(ligneSplitted[13])+"]"])
                localPDF.extend(["sigmaF["+str(ligneSplitted[14])+"]"])
                localPDF.extend(["vFrac["+str(ligneSplitted[15])+"]"])
                localPDF.extend(["wf["+str(ligneSplitted[16])+"]"])
                localPDF.extend(["wp["+str(ligneSplitted[17])+"]"])        
                localPDF.extend(TemplateSignal_Landau)
                nomPDF=ID+"_pt_"+ligneSplitted[0]+"_"+ligneSplitted[1]
                setattr(process.TnP_ElecID.PDFs, nomPDF,localPDF)
                continue
            #if (("absSCeta_bin2" in ligneSplitted[0]) or ("absSCeta_bin1" in ligneSplitted[0])):
            if (("absSCeta_bin2" in ligneSplitted[0])):
                localPDF = cms.vstring()
                localPDF.extend(["lef["+str(ligneSplitted[5])+"]"])
                localPDF.extend(["lep["+str(ligneSplitted[6])+"]"])
                localPDF.extend(["mean1f["+str(ligneSplitted[7])+"]"])
                localPDF.extend(["mean1p["+str(ligneSplitted[8])+"]"])
                localPDF.extend(["mean2f["+str(ligneSplitted[9])+"]"])
                localPDF.extend(["mean2p["+str(ligneSplitted[10])+"]"])
                localPDF.extend(["sigma1f["+str(ligneSplitted[12])+"]"])
                localPDF.extend(["sigma1p["+str(ligneSplitted[13])+"]"])
                localPDF.extend(["sigma2f["+str(ligneSplitted[14])+"]"])
                localPDF.extend(["sigma2p["+str(ligneSplitted[15])+"]"])
                localPDF.extend(["vFrac["+str(ligneSplitted[16])+"]"])
                localPDF.extend(["vPropF["+str(ligneSplitted[17])+"]"])
                localPDF.extend(["vPropP["+str(ligneSplitted[18])+"]"])
                localPDF.extend(TemplateSignal_BadZone)
                nomPDF=ID+"_pt_"+ligneSplitted[0]+"_"+ligneSplitted[1]
                setattr(process.TnP_ElecID.PDFs, nomPDF,localPDF)
                continue
            localPDF = cms.vstring()
            localPDF.extend(["alpha["+str(ligneSplitted[2])+"]"])
            localPDF.extend(["alphaF["+str(ligneSplitted[3])+"]"])
            localPDF.extend(["lef["+str(ligneSplitted[7])+"]"])
            localPDF.extend(["lep["+str(ligneSplitted[8])+"]"])
            localPDF.extend(["mean["+str(ligneSplitted[9])+"]"])
            localPDF.extend(["n["+str(ligneSplitted[10])+"]"])
            localPDF.extend(["scale["+str(ligneSplitted[12])+"]"])
            localPDF.extend(["sigma["+str(ligneSplitted[13])+"]"])
            localPDF.extend(["sigmaF["+str(ligneSplitted[14])+"]"])
            localPDF.extend(["vFrac["+str(ligneSplitted[15])+"]"])
            localPDF.extend(TemplateSignal)
            nomPDF=ID+"_pt_"+ligneSplitted[0]+"_"+ligneSplitted[1]
    # print nomPDF
    #   print localPDF
            setattr(process.TnP_ElecID.PDFs, nomPDF,localPDF)

print process.TnP_ElecID.PDFs



ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_FO" or "_Preselec" or "_IP" or "_Convs" or "_FOnoIso" in scenario: ALLBINS = [("ptFO",PT_ETA_BINS_POG)]
if "_NM1ISO" in scenario: ALLBINS = [("ptNM1ISO",PT_ETA_BINS_POG_NM1ISO)]
if "_NM1ID" in scenario: ALLBINS = [("ptNM1ID",PT_ETA_BINS_POG_NM1ID)]
if "_NM1FO" in scenario: ALLBINS = [("ptNM1FO",PT_ETA_BINS_POG_NM1FO)]
if "_NM1DETISO" in scenario: ALLBINS = [("ptNM1DETISO",PT_ETA_BINS_POG_NM1DETISO)]
if "_isoOnly" in scenario: ALLBINS = [("ptIsoOnly",PT_ETA_BINS_POG_ISOONLY)]
if "_PSEL" in scenario: ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_NM1IP" in scenario: ALLBINS = [("ptNM1IP",PT_ETA_BINS_POG_NM1IP)]
if "_NM1CONVS" in scenario: ALLBINS = [("ptNM1CONVS",PT_ETA_BINS_POG_NM1CONV)]
if "_NM1PSEL" in scenario: ALLBINS = [("ptNM1PSEL",PT_ETA_BINS_POG_NM1PRESEL)]
if "_ALLNOPSEL" in scenario: ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_ALLNOCONV" in scenario: ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_ALLNOIP" in scenario: ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_ALLTHESELEC" in scenario: ALLBINS = [("pt",PT_ETA_BINS_POG)]
if "_ALLTHESELEC2" in scenario: ALLBINS = [("pt2",PT_ETA_BINS_POGL1)]
if "_ALLTHESELEC3" in scenario: ALLBINS = [("pt3",PT_ETA_BINS_POGL2)]
if "_ALLTHESELECTRIG" in scenario: ALLBINS = [("ptSelecTrig",PT_ETA_BINS_POG_TRIG)]
if "_ALLTHESELECTRIG2" in scenario: ALLBINS = [("ptSelecTrig2",PT_ETA_BINS_POG_TRIG2)]
if "_ALLTHESELECTRIG3" in scenario: ALLBINS = [("ptSelecTrig3",PT_ETA_BINS_POG_TRIG3)]



if len(args) > 1 and args[1] not in IDS: IDS += [ args[1] ]
for ID in IDS:
    if len(args) > 1 and ID != args[1]: continue
    if "forBg" in scenario:
        ALLBINS = [("pt",PT_ETA_BINS_POG_BG)]
    print ALLBINS
    # ALLBINS += [("ptEta_smurfs", PT_ETA_BINS_SMURFS)]
    #if "passElec_FO_ISO_ID" in ID: ALLBINS += [("vtx",VTX_BINS)]
    #if "Pog" in ID: ALLBINS += [("ptPOG",PT_ETA_BINS_POG)]
    for X,B in ALLBINS:
        print "coucou on va faire ", ID , " avec les bins ", X
        if len(args) > 2 and X not in args[2:]: continue
        module = process.TnP_ElecID.clone(OutputFileName = cms.string("TnP_ElecID_%s_%s_%s.root" % (scenario, ID, X)))
        #shape = cms.vstring("CBBreitWignerPlusExponentialBackground","*pt_bin5*","CBBreitWignerPlusPolyBackground")
        #   shape = cms.vstring("SignalFromMC","*pt_bin5*","SignalFromMC_Landau")
        shape = cms.vstring("CBBreitWignerPlusExponentialBackground")
        if "mc" in scenario:
            shape = cms.vstring("SignalFromMC","*absSCeta_bin0*pt_bin5*","SignalFromMC_Landau","*absSCeta_bin1*pt_bin5*","SignalFromMC_Landau","*absSCeta_bin1*pt_bin0*","SignalFromMC","*absSCeta_bin1*pt_bin1*","SignalFromMC","*absSCeta_bin1*pt_bin2*","SignalFromMC","*absSCeta_bin1*pt_bin3*","SignalFromMC","*absSCeta_bin1*pt_bin4*","SignalFromMC","*absSCeta_bin2*","SignalFromMC_BadZone","*absSCeta_bin3*pt_bin5*","SignalFromMC_Landau","*absSCeta_bin4*pt_bin5*","SignalFromMC_Landau")
        if "data" in scenario:
            for i in range(5):
                for j in range(6):
                    shape.extend(["*absSCeta_bin"+str(i)+"*pt_bin"+str(j)+"*",ID+"_pt_absSCeta_bin"+str(i)+"_pt_bin"+str(j)])
        if "forBg" in scenario:
            shape = cms.vstring("SignalFromMCMyBgShape")
        print shape
        DEN = B.clone(); num = ID;
        if "_from_" in ID:
            parts = ID.split("_from_")
            num = parts[0]
            setattr(DEN, parts[1], cms.vstring("pass"))
        if scenario.find("tagiso") != -1:  
            DEN.tag_combRelIso = cms.vdouble(-1, 0.1)
        if scenario.find("loosetagiso") != -1:  
            DEN.tag_combRelIso = cms.vdouble(-1, 0.2)
        if scenario.find("probeiso") != -1:
            DEN.isoTrk03Abs = cms.vdouble(-1, 3)
        #if scenario.find("calo") != -1: DEN.caloCompatibility = cms.vdouble(0.9,1.1)  # same as above, I think.
        if "mc" in scenario:
            if num == "Mu24": num = "Mu15"
            if num == "IsoMu17": num = "IsoMu15"
            if num == "DoubleMu7": num = "DoubleMu3"
            if num == "Mu8_forEMu": num = "DoubleMu3"
            if num == "Mu17_forEMu": num = "DoubleMu3"
        if "EG5" in scenario: DEN.pair_nL1EG5 = cms.vdouble(0.5,999)
        
        if "data" in scenario:
            setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
                                                            EfficiencyCategoryAndState = cms.vstring(num,"pass"),
                                                            UnbinnedVariables = cms.vstring("mass"),
                                                            BinnedVariables = DEN,
                                                            BinToPDFmap = shape#cms.vstring(shape)
                                                            ))
        if "RunA" in scenario: theWeight="weight_runA"
        elif "RunB" in scenario: theWeight="weight_runB"
        elif "RunC" in scenario: theWeight="weight_runC"
        elif "RunD" in scenario: theWeight="weight_runD"
        if scenario.find("mc") != -1:
            setattr(module.Efficiencies, ID+"_"+X, cms.PSet(
                                                            EfficiencyCategoryAndState = cms.vstring(num,"pass"),
                                                            UnbinnedVariables = cms.vstring("mass",theWeight),
                                                            BinnedVariables = DEN,
                                                            BinToPDFmap = shape
                                                            ))
        setattr(process, "TnP_ElecID_"+ID+"_"+X, module)        
        setattr(process, "run_"+ID+"_"+X, cms.Path(module))

