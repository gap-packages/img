---
layout: default
---

# GAP Package {{site.data.package.name}}

{{site.data.package.abstract}}

The current version of this package is version {{site.data.package.version}}.
For more information, please refer to [the package manual]({{site.data.package.doc-html}}).
There is also a [README](README.html) file.

## Dependencies

This package requires GAP version {{site.data.package.GAP}}
{% if site.data.package.needed-pkgs %}
The following other GAP packages are needed:
{% for pkg in site.data.package.needed-pkgs %}
- {% if pkg.url %}<a href="{{ pkg.url }}">{{ pkg.name }}</a>{% else %}{{ pkg.name }}{% endif %} {{ pkg.version }}{% endfor %}
{% endif %}
{% if site.data.package.suggested-pkgs %}
The following additional GAP packages are not required, but suggested:
{% for pkg in site.data.package.suggested-pkgs %}
- {% if pkg.url %}<a href="{{ pkg.url }}">{{ pkg.name }}</a>{% else %}{{ pkg.name }}{% endif %} {{ pkg.version }}{% endfor %}
{% endif %}

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
