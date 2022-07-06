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
            for (let i in lines) {
                if (i == 0) {continue;}
                if (/^    [\w_. -]+: /.test(lines[i])) {
                    let br = (i == 1) ? "" : "<br />";
                    items.push(br + lines[i].replace(
                        /^    ([\w_. -]+): /,
                        '<code class="mkapi-item-name">$1</code>'
                        + '<span class="mkapi-item-dash"> — </span>'
                    ));
                } else if (/^\s*-.+/.test(lines[i])) {
                    items.push("<br />" + lines[i]);
                } else {
                    items.push(lines[i]);
                }
            }
            p.addClass("proc-attr").html(
                '<div class="mkapi-section">'
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
