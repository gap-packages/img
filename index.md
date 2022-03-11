---
layout: default
---

# GAP Package {{site.data.package.name}}

{{site.data.package.abstract}}.

The current version of this package is version {{site.data.package.version}}.
For more information, please refer to [the package manual]({{site.data.package.doc-html}}).

## Dependencies

This package requires GAP version {{site.data.package.GAP}}
{% if site.data.package.needed-pkgs %}
The following GAP packages are needed:
{% for pkg in site.data.package.needed-pkgs %}
- [FR](https://gap-packages.github.io/fr), version at least 2.4.0.

The following additional GAP packages are not required, but suggested:
- [Float](https://gap-packages.github.io/float/), version at least 0.9.0. It is only necessary for high-precision (greater than 64-bit IEEE754) calculations

Note also that IMG interacts with a GUI that displays interactively Julia sets on the sphere; it is included in the tarball, but has to be added manually if you clone the git IMG repository. Do this with `git clone https://github.com/laurentbartholdi/rsserver.git` in the root directory of IMG.

## Installation (MacOS)

Installation has been tested with [homebrew](https://brew.sh/) and [macports](https://www.macports.org/). You definitely will need at least one of them. It is assumed that you already have a running installation of GAP, either packaged or self-compiled.

* Download or clone the latest version of this package in $HOME/Library/Preferences/GAP/pkg (you may need to create the directory beforehand)
* For homebrew: `brew install cmake suite-sparse imagemagick node`; for macports: `port install cmake SuiteSparse ImageMagick nodejs17`
* go to the root of the package ($HOME/Library/Preferences/GAP/pkg/img for example)
* `./configure --with-gaproot=<root> --with-cholmod=<cholmod-root>`; here `<root>` is typically `/opt/local/Cellar/gap/4.11.1/libexec` and `<cholmod-root>` is typically `/opt/local`
* `make`

## Installation (Linux)

(to be inserted from Lukas Geyer's notes)

## Author{% if site.data.package.authors.size != 1 %}s{% endif %}
{% for person in site.data.package.authors %}
{% if person.url %}<a href="{{ person.url }}">{{ person.name }}</a>{% else %}{{ person.name }}{% endif %}{% unless forloop.last %}, {% endunless %}{% else %}
{% endfor %}

{% if site.github.issues_url %}
## Feedback

For bug reports, feature requests and suggestions, please use the
[issue tracker]({{site.github.issues_url}}).
{% endif %}
