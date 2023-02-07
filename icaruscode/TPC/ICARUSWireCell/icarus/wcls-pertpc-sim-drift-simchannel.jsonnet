// Same configuration as in wcls-sim-drift-simchannel.jsonnet
// except that this produces four instances of std::vector<RawDigits>
// one per physics module (WW, WE, EE, EW) in ICARUS

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
// local params = import 'pgrapher/experiment/icarus/simparams.jsonnet';
local base = import 'pgrapher/experiment/icarus/simparams.jsonnet';
local params = base {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    // drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
  files: super.files {
    fields: [ std.extVar('files_fields'), std.extVar('files_fields'), std.extVar('files_fields'), std.extVar('files_fields'),
              std.extVar('files_fields'), std.extVar('files_fields'), std.extVar('files_fields'), std.extVar('files_fields')],
  },
};
local volname = ["EE", "EW", "WE", "WW"];
local cryo_number = std.extVar('cryo_number');

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/icarus/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local output = 'wct-sim-ideal-sig.npz';


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  // depos: wcls.input.depos(name="", art_tag="ionization"),
  depos: wcls.input.depos(name='electron', art_tag='ionization'),  // default art_tag="blopper"
};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

// A ``duo'' anode consists of two ``splits''
local duoanodes = {
  type: 'MegaAnodePlane',
  name: 'duoanode%d' %cryo_number,
  data: {
    anodes_tn: [wc.tn(a) for a in tools.anodes[2*cryo_number:2*(cryo_number+1)]],
  }, 
};
local wcls_output = {
  // ADC output from simulation
  // sim_digits: wcls.output.digits(name="simdigits", tags=["orig"]),

//  sim_digits: g.pnode({
//    type: 'wclsFrameSaver',
//    name: 'simdigits%d' %cryo_number,
//    data: {
//      // anode: wc.tn(tools.anode),
//      // anode: wc.tn(duoanodes),
//      anode: wc.tn(mega_anode),
//      digitize: true,  // true means save as RawDigit, else recob::Wire
//      //frame_tags: ['daq%d' %cryo_number],
//      frame_tags: ['TPC%s' %volname[cryo_number]],
//      // Three options for nticks:
//      // - If nonzero, force number of ticks in output waveforms.
//      // - If zero, use whatever input data has. (default)
//      // - If -1, use value as per LS's detector properties service.
//      // nticks: params.daq.nticks,
//      // nticks: -1,
//      // chanmaskmaps: ['bad'],
//    },
//  // }, nin=1, nout=1, uses=[duoanodes]),
//  }, nin=1, nout=1, uses=[mega_anode]),

  sim_digits: [ 
  g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits%d' %n,
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(duoanodes),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      //frame_tags: ['daq%d' %n],
      frame_tags: ['TPC%s' %volname[n]],
      // Three options for nticks:
      // - If nonzero, force number of ticks in output waveforms.
      // - If zero, use whatever input data has. (default)
      // - If -1, use value as per LS's detector properties service.
      // nticks: params.daq.nticks,
      // nticks: -1,
      // chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[duoanodes])
  for n in std.range(0,3)],

  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: wcls.output.digits(name='nfdigits', tags=['raw']),

  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),

  // save "threshold" from normal decon for each channel noise
  // used in imaging
  sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
// local sn_pipes = sim.splusn_pipelines;
local analog_pipes = sim.analog_pipelines;

local perfect = import 'pgrapher/experiment/icarus/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];


// local nf_maker = import 'pgrapher/experiment/icarus/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/icarus/sp.jsonnet';
local sp = sp_maker(params, tools);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
  name: 'postdrift%d' %cryo_number,
  data: {
    artlabel: 'simpleSC%s' %volname[cryo_number],  // where to save in art::Event
    // anodes_tn: [wc.tn(duoanodes)],
    // anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    anodes_tn: [wc.tn(a) for a in tools.anodes[2*cryo_number:2*(cryo_number+1)]],
    rng: wc.tn(rng),
    tick: params.daq.tick,
    start_time: -0.34 * wc.ms, // TriggerOffsetTPC from detectorclocks_icarus.fcl
    readout_time: params.daq.readout_time,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: 100 * wc.mm,
    v_to_rp: 100 * wc.mm,
    y_to_rp: 100 * wc.mm,

    // GP: The shaping time of the electronics response (2.2us) shifts the peak
    //     of the field response time. Eyeballing simulation times, it does this
    //     by a bit less than the 2.2us.
    //
    //     N.B. for future: there is likely an additional offset on the two induction
    //     planes due to where the deconvolution precisely defines where the "peak"
    //     of the pulse is. One may want to refine these parameters to account for that.
    //     This perturbation shouldn't be more than a tick or two.
    u_time_offset: 2.0 * wc.us,
    v_time_offset: 2.0 * wc.us,
    y_time_offset: 2.0 * wc.us,

    g4_ref_time: -1500 * wc.us, // G4RefTime from detectorclocks_icarus.fcl
    use_energy: true,
  },
// }, nin=1, nout=1, uses=[duoanodes]);
}, nin=1, nout=1, uses=tools.anodes);
// }, nin=1, nout=1, uses=tools.anodes[2*cryo_number:2*(cryo_number+1)]);

local nicks = ["incoTPCEE","incoTPCEW","incoTPCWE","incoTPCWW", "coheTPCEE","coheTPCEW","coheTPCWE","coheTPCWW"];
local scale_int = std.extVar('int_noise_scale');
local scale_coh = std.extVar('coh_noise_scale');
local model_int = {
        type: "GroupNoiseModel",
        name: nicks[cryo_number],
        data: {
            // This can also be given as a JSON/Jsonnet file
            spectra: params.files.noisegroups[cryo_number],
            groups: params.files.wiregroups,
            scale: scale_int,
            nsamples: params.daq.nticks,
            tick: params.daq.tick,
        }
    };
local model_coh = {
        type: "GroupNoiseModel",
        name: nicks[4+cryo_number],
        data: {
            // This can also be given as a JSON/Jsonnet file
            spectra: params.files.noisegroups[4+cryo_number],
            groups: params.files.wiregroups,
            scale: scale_coh,
            nsamples: params.daq.nticks,
            tick: params.daq.tick,
        }
    };



local add_noise = function(model, n,t) g.pnode({
    type: t,
    name: "addnoise%d-" %n + model.name,
    data: {
        rng: wc.tn(tools.random),
        dft: wc.tn(tools.dft),
        model: wc.tn(model),
        nsamples: params.daq.nticks,
    }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]);
local noises = add_noise(model_int, cryo_number,"IncoherentAddNoise");
local coh_noises = add_noise(model_coh, 4 + cryo_number,"CoherentAddNoise");

// local digitizer = sim.digitizer(mega_anode, name="digitizer", tag="orig");
local digitizers = sim.digitizer(mega_anode, name="digitizer%d-" %cryo_number + mega_anode.name, tag="TPC%s"%volname[cryo_number]);
// local digitizers = [
//     sim.digitizer(mega_anode, name="digitizer%d-" %n + mega_anode.name, tag="TPC%s"%volname[n])
//     for n in std.range(0,3)];

local frame_summer = g.pnode({
        type: 'FrameSummer',
        name: 'framesummer%d' %cryo_number,
        data: {
            align: true,
            offset: 0.0*wc.s,
        },
    }, nin=2, nout=1);

local actpipe = g.pipeline([noises, coh_noises, digitizers, wcls_output.sim_digits[cryo_number]], name="noise-digitizer%d" %cryo_number);
// local actpipe = g.pipeline([noises, coh_noises, digitizers[cryo_number], wcls_output.sim_digits], name="noise-digitizer%d" %cryo_number);
local util = import 'pgrapher/experiment/icarus/funcs.jsonnet';
local outtags = ['orig%d' % cryo_number];
local pipe_reducer = util.minifansummer('DepoSetFanout', analog_pipes[2*cryo_number:2*(cryo_number+1)], frame_summer, actpipe, cryo_number, 'FrameFanin', 'minifansummer', outtags);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

// local graph = g.pipeline([wcls_input.depos, drifter,  wcls_simchannel_sink, bagger, pipe_reducer, retagger, wcls_output.sim_digits, sink]);
local graph = g.pipeline([wcls_input.depos, drifter,  wcls_simchannel_sink, bagger, pipe_reducer, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
