// Hide base attributes for procs

$('div[id^="biopipen.ns."]').each(function() {
    let node = $(this);
    let parts = node.attr('id').split('.');
    if (parts.length == 3) {
        return;
    }
    node.find(".mkapi-section.attributes").remove();
    node.find(".mkapi-section.classes").remove();
    node.find(".mkapi-section.methods").remove();
    node.find(".mkapi-members").remove();

    // format items
    node.find(".mkapi-section-body > p").each(function() {
        let p = $(this);
        let lines = p.html().split("\n");
        if (/^(Input|Output|Envs|Requires):$/.test(lines[0])) {
            let items = [];
            console.log(lines);
            for (let i = 0; i < lines.length; i++) {
                if (i === 0) { continue; }
                if (lines[0] === "Requires:" && /^\s+-? \w+:\s*(\| )?/.test(lines[i])) {
                    let key = lines[i].replace(/^\s+(-? \w+:\s*(?:\| )?).+/, "$1");
                    let val = lines[i].replace(/^\s+-? \w+:\s*(?:\| )?/, "");
                    // console.log(key, "-", val);
                    let kls = key[0] === "-" ? " mkapi-item-first" : "";
                    items.push(
                        '<code class="mkapi-item-name' + kls + '">'
                        + key.replace(/(^[\s-]+|[:\|\s]+$)/g, "")
                        + '</code>'
                        + '<span class="mkapi-item-dash"> — </span>'
                        + '<code>' + val + '</code>'
                        + '<br />'
                    );
                } else if (/^    [\w_. -]+: /.test(lines[i])) {
                    let br = (i === 1) ? "" : "<br />";
                    items.push(br + lines[i].replace(
                        /^    ([\w_.-]+): /,
                        '<code class="mkapi-item-name">$1</code>'
                        + '<span class="mkapi-item-dash"> — </span>'
                    ));
                } else if (/^\s*-.+/.test(lines[i])) {
                    items.push("<br />" + lines[i]);
                } else if (lines[0] === "Requires:") {
                    items.push('<code class="pipen-requires-check-code">' + lines[i] + "</code>");
                } else {
                    items.push(lines[i]);
                }
            }
            let require_kls = (lines[0] === "Requires:") ? " pipen-requires" : "";
            p.addClass("proc-attr").html(
                '<div class="mkapi-section' + require_kls + '">'
                + '<div class="mkapi-section-name">'
                + '<span class="mkapi-section-name-body">'
                + lines[0]
                + '</span>'
                + '</div>'
                + '<div class="mkapi-section-body">'
                + items.join("\n")
                + '</div>'
                + '</div>'
            );
        }
    });
});

// TOC
$('div.md-sidebar--secondary nav.md-nav--secondary nav.md-nav a.md-nav__link[href^="#biopipenns"]')
    .next()
    .remove();
