#ifndef KALOS_TRAMO_OUTLIER_ROUNDS_HPP
#define KALOS_TRAMO_OUTLIER_ROUNDS_HPP

#include "kalos/tramo_outlier_scan.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace kalos {

// Item 23-F F.2: the FAITHFUL multi-round TRAMO outlier state machine, the
// FastOutliersDetector.process port of the pinned jdplus 3.9.0 stack (the
// round-0 decision unit of item 23-B observes the first scan only). The
// driver shares its kernels with that unit (src/tramo_scan_machine.hpp); the
// bytecode pin and the per-round corpus (37 cases of full field dumps,
// byte-identical dual-machine captures) live in
// docs/contracts/tramo-outlier-scan-contract.md and
// validation/evidence/native/tramo-outlier-rounds/.
//
// The result trace mirrors the probe's dump sequence 1:1 (init, then one
// event per estimateModel / loop-head verifyModel / scan, the add and cap
// records, the conditional tail estimateModel, the drain loop's verify /
// estimate pairs and the final manret / manouts records), so a corpus case
// replays field-by-field. All doubles are bit-exact against the witness
// captures.
//
// Caps use the FastOutliersDetector FIELD semantics. The pinned
// OutliersDetectionModule builder hands its maxOutliers and maxRound values
// through a constructor whose parameter order lands the FIRST builder value
// on field maxRound and the SECOND on field maxOutliers (the 23-F F.2 javap
// decode); the trace and this unit always speak the FIELD values, so the
// pinned defaults are (max_round, max_outliers) = (30, 50) either way and a
// corpus case with custom caps records the fields already crossed in its
// witness line.
//
// Request validation is the 23-B gate verbatim (prefix tramo_outlier_rounds)
// and the machine-internal dead ends (a rejected size gate, "Ols failed",
// "Gls failed", the scan table without admissible cells, the JVM exception
// paths of the null-field reads and of a failing driver read) map onto
// Status::Failed with the driver's exit note, never onto an approximation.

// One probe state dump (FastOutliersDetector fields + driver views). The
// *_null flags carry the JVM null fields; the witness prints them as null,
// vectors as %.17g#hex bit patterns.
struct TramoOutlierRoundState {
    int round;
    bool rflag;             // residual-flag (OLS branch)
    bool estimation_step;   // estimationStep
    bool backward_step;     // backwardStep
    bool exit;              // exit (set on re-removing lastremoved)
    std::vector<std::array<int, 2>> outliers; // {position, kind} in append order
    bool lastremoved;       // the witness prints the pair's length (2) or null
    std::vector<double> coeff;
    bool coeff_null;
    std::vector<double> tstats;
    bool tstats_null;
    std::vector<double> res;
    bool res_null;
    std::vector<double> pars; // ARIMA parameters at JVM coordinates
    std::size_t xcount;       // user columns of the driver's RegArimaModel
};

// One SingleOutlierDetector probe dump (dumpSod fields). allowed_counts is
// the per-column count of live cells of the TableOfBoolean (the mask
// persistence pin); maxt_threw mirrors the searchMax (-1, -1) table read.
struct TramoOutlierRoundScanDump {
    bool el_null;
    std::size_t el_count;
    double mad;
    int lower_bound, upper_bound; // the detector's effective bounds
    bool t_table_present;
    bool allowed_present;
    std::size_t allowed_rows, allowed_cols;
    std::vector<std::size_t> allowed_counts;
    bool maxt_threw;
    double max_t;
    bool pick;                    // a winner exists (searchMax fired)
    int pick_position, pick_type;
    double coefficient;           // the coefficient table cell of the winner
};

// One trace record. Tags follow the witness names: "init", "est", "verify",
// "scan" (state + sod), "add" (state + add_* fields), "tailest", "drainv",
// "draine", "drainend", "manret" and "manouts". All outlier "types" (picks,
// add_*, state/outliers) are the searchMax oMax, the MIX COLUMN index - the
// kind code only in the full AO,LS,TC,SO mix.
struct TramoOutlierRoundEvent {
    std::string tag;
    bool has_state;
    TramoOutlierRoundState state;
    bool has_sod;
    TramoOutlierRoundScanDump sod;
    // The add line of the witness (valid with tag "add").
    int add_position, add_type;
    double add_coefficient;
    // The manret line (valid with tag "manret").
    bool manret_returned;
    std::string note; // run | est | exit | sodfalse | gate | cap | roundcap
                      // | maxtthrew | stepcap | draincap
    int manret_round;
};

struct TramoOutlierRoundsResult {
    // Ok: process returned true (the selection is final). Failed: process
    // returned false - the module reports ProcessingResult.Failed.
    enum class Status { Ok, Failed };
    Status status;
    std::string note;  // the exit note of the final manret record
    int round;         // the round field at exit
    // The final outlier set (driver getOutliers; empty on most Failed runs).
    std::vector<std::array<int, 2>> outliers;
    std::vector<TramoOutlierRoundEvent> trace;
};

// Stationary regular AR(1)/AR(2) and seasonal ARMA(1,1) paired with
// regular AR(2) or ARMA(1,1), with a mean, a complete response
// and the AO/LS/TC factory family qualified against the JVM round trace.
bool tramo_stationary_scan_supported(const TramoOutlierScanRequest& request);

// Runs the full multi-round machine on the 23-B request shape (the caps are
// read as FIELD values, see the header note; the scan unit ignores them).
// Static request violations throw std::invalid_argument exactly like
// scan_tramo_outliers; legal requests never throw for pinned numerical dead
// ends, they report Status::Failed.
TramoOutlierRoundsResult tramo_outlier_rounds(
    const TramoOutlierScanRequest& request);

} // namespace kalos
#endif
