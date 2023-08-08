use clap::arg;
use clap::Command;
use jam_rs::{compare, file_io};

fn cli() -> Command {
    Command::new("jam")
        .about("Just another minhasher, obviously blazingly fast")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .allow_external_subcommands(true)
        .subcommand(
            Command::new("sketch")
                .about("Create a sketch for file(s)")
                .arg(arg!(-i --input <FILE> "Input file"))
                .arg_required_else_help(true)
                .arg(arg!(-o --output <FILE> "Output file"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("diff")
                .about("Compare two commits")
                .arg(arg!(base: [COMMIT]))
                .arg(arg!(head: [COMMIT]))
                .arg(arg!(path: [PATH]).last(true))
                .arg(
                    arg!(--color <WHEN>)
                        .value_parser(["always", "auto", "never"])
                        .num_args(0..=1)
                        .require_equals(true)
                        .default_value("auto")
                        .default_missing_value("always"),
                ),
        )
        .arg(arg!(--version).exclusive(true))
        .arg(arg!(--threads -t <THREADS>).global(true))
}

fn main() {
    let matches = cli().get_matches();

    match matches.subcommand() {
        Some(a) => (),
        None => (),
    }
}
