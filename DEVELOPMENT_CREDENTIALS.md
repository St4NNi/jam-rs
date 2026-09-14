# Limited-credentials development build

This checkout's selected default development build restricts S3 credentials. It is
an explicit development capability change, not a release compatibility promise.
Local indexed tracing is unchanged. No cloud deployment or release is implied.

Use complete `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` environment values,
with optional `AWS_SESSION_TOKEN` or `AWS_SECURITY_TOKEN`, or static credentials
from `AWS_SHARED_CREDENTIALS_FILE` (otherwise `~/.aws/credentials`). `AWS_PROFILE`,
then `AWS_DEFAULT_PROFILE`, selects a section; otherwise `[default]` is used.
Static profile session/security tokens are preserved. A partial environment
selection fails instead of trying another profile or account. No anonymous
fallback is used. These static credentials do not refresh automatically.

STS web identity, container credential endpoints and EC2 instance metadata are
unavailable. Explicit web identity, role, container or metadata endpoint environment
settings fail with `unsupported credential provider` before constructing the S3
configuration or making a request, even when static credentials are also set.
Remove such settings only when intentionally selecting static credentials.
Role-only profiles and missing/invalid static credentials return the same clear
error. Before selecting profile keys, both the selected shared-credentials section
and its `AWS_CONFIG_FILE` (otherwise `~/.aws/config`) section are checked for role,
source-profile, credential-source, web-identity, credential-process and `sso_*`
directives. These fail even when static base-account keys are present; silently
using those keys could select the wrong account. Nonselected profiles do not
block valid selected static credentials. Explicit complete environment credentials
retain precedence over profile files. AWS config-file role chains, credential
processes and SSO are unsupported. The existing locked `rust-ini` and `home`
packages are direct dependencies so this check uses the upstream INI syntax and
home-path behavior, with no second parser or new dependency package.

The manifest selects `rust-s3 0.37.2` with `fail-on-err,sync`, and selects the
existing `attohttpc 0.30` transport's `tls-rustls` feature directly, with defaults
disabled on both edges. It does not select `rust-s3/sync-rustls-tls`, because that
also selects `aws-creds/rustls-tls` and hence `aws-creds/http-credentials`.
TLS certificate verification, AWS request signing, session-token headers, HEAD
metadata checks, ranged GET validation, and error redaction remain enabled.

Build with the checked-in manifest and lockfile using `cargo build --release
--locked` and the recorded external target/temp roots and generic x86-64 flags.
Verify `cargo tree --locked -e features -i aws-creds` has no enabled HTTP-provider
or provider-TLS features, and `cargo tree --locked -e features -i attohttpc` retains
transport TLS. Cargo features are additive: another consumer can re-enable an
upstream feature, so a library consumer must verify its complete resolved graph.
There are no local patches to transfer. A future registry installation must use a
published package containing this manifest, its default feature selection and its
lockfile (`cargo install jam-rs --version <that-version> --locked`). Installing an
older published version does not reproduce this profile; publication is outside
this development checkpoint.

## Dependency findings and scoped review

Both upstream packages unconditionally depend on `quick-xml ^0.38`; the locked
version is 0.38.4. Official upstream manifests inspected on 2026-09-14 still select
this series. The patched version for both findings below is at least 0.41.0,
outside those compatible dependency requirements. No vendor patch, blanket ignore,
or unrelated dependency update is used.

| Finding | Development disposition |
| --- | --- |
| [RUSTSEC-2026-0194](https://rustsec.org/advisories/RUSTSEC-2026-0194.html), duplicate-attribute quadratic CPU use | Remains a dependency finding; not patched or ignored. |
| [RUSTSEC-2026-0195](https://rustsec.org/advisories/RUSTSEC-2026-0195.html), namespace allocation denial of service | Remains a dependency finding; not patched or ignored. |
| [RUSTSEC-2024-0436](https://rustsec.org/advisories/RUSTSEC-2024-0436.html), unmaintained `paste 1.0.15` | Narrow informational exception for this development build only. |

Jam's enabled S3 call sites are `range_source::S3Config::bucket`,
`Bucket::head_object`, and `Bucket::get_object_range`. HEAD converts headers to
metadata; ranged GET returns bytes. The synchronous `fail-on-err` backend converts
HTTP error bodies to text, and Jam exposes its fixed error message. None of these
paths calls the XML deserializer. The upstream STS XML path is inside the removed
`http-credentials` feature. Other XML entry points remain in rust-s3 (including
listing, attributes, location, CORS, lifecycle, and multipart operations), but Jam
does not call them. This is a scoped call-path audit, not proof that the package
has disappeared or an advisory-clean claim. Local responders exercise malformed
XML-shaped successful data and error responses, signed requests and session
tokens; tests cannot establish general XML parser safety.

The paste exception applies only to `jam-rs -> statrs 0.18.0 default/nalgebra ->
nalgebra 0.33.2 std -> simba 0.9.1 -> paste 1.0.15 default` (also reached through
statrs' default rand features). Paste is a build-time procedural macro with an
unmaintained informational advisory and no patched version. Replacing the math
dependency chain is unrelated to this bounded tracing change, so that maintenance
risk is accepted for this development build. Revisit at the next dependency review
and before release, or earlier if the path changes or a vulnerability or supported
replacement/update is identified. This exception does not cover either quick-xml
denial-of-service finding. No audit ignore is added.

Upstream contracts: [rust-s3 manifest](https://raw.githubusercontent.com/durch/rust-s3/v0.37.2/s3/Cargo.toml),
[aws-creds manifest](https://raw.githubusercontent.com/durch/rust-s3/v0.37.2/aws-creds/Cargo.toml),
[Cargo feature unification](https://doc.rust-lang.org/cargo/reference/resolver.html#features).
