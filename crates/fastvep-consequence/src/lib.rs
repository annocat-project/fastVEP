mod predictor;
mod splice;
pub mod sv_predictor;

pub use predictor::{
    vep_input_position, AlleleConsequenceResult, ConsequencePredictor, PredictionResult, TranscriptConsequence,
};
